#include <iostream>
#include <vector>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <sndfile.h>
#include <fftw3.h>
#include <Eigen/Dense>

// ========== Config ==========
constexpr int   SAMPLE_RATE   = 16000;
constexpr int   FRAME_LENGTH  = 400;     // 25 ms
constexpr int   HOP_LENGTH    = 160;     // 10 ms
constexpr int   N_FFT         = 512;     // padded FFT
constexpr int   N_MELS        = 128;
constexpr int   TARGET_FRAMES = 1024;

constexpr float F_MIN         = 20.0f;
constexpr float F_MAX         = SAMPLE_RATE / 2.0f;
constexpr float PREEMPHASIS   = 0.97f;
constexpr float MEL_FLOOR     = 1.192092955078125e-07f;  // same as HF path
constexpr float MEAN          = -4.2677393f;
constexpr float STD           = 4.5689974f;

// ========== Kaldi mel scale ==========
inline float hz_to_mel_kaldi(float hz) {
    return 1127.0f * std::log(1.0f + hz / 700.0f);
}
inline float mel_to_hz_kaldi(float mel) {
    return 700.0f * (std::exp(mel / 1127.0f) - 1.0f);
}

// Build mel filterbank like Kaldi / torchaudio.functional.get_mel_banks
// Triangularization happens in MEL space

Eigen::MatrixXf build_mel_filterbank() {
    const int n_fft_bins = N_FFT / 2 + 1;
    Eigen::MatrixXf fb = Eigen::MatrixXf::Zero(N_MELS, n_fft_bins);

    const float mel_min = hz_to_mel_kaldi(F_MIN);
    const float mel_max = hz_to_mel_kaldi(F_MAX);
    const float mel_step = (mel_max - mel_min) / (N_MELS + 1);

    // FFT bin center freqs (Hz) → mel
    const float fft_bin_width = float(SAMPLE_RATE) / float(N_FFT);
    std::vector<float> mel_of_bin(n_fft_bins);
    for (int k = 0; k < n_fft_bins; ++k) {
        float hz = k * fft_bin_width;
        mel_of_bin[k] = hz_to_mel_kaldi(hz);
    }

    // For each mel filter, define left/center/right in mel,
    // then compute triangular weights against mel_of_bin
    for (int m = 0; m < N_MELS; ++m) {
        float left_m   = mel_min +  m      * mel_step;
        float center_m = mel_min + (m + 1) * mel_step;
        float right_m  = mel_min + (m + 2) * mel_step;

        for (int k = 0; k < n_fft_bins; ++k) {
            float mel_k = mel_of_bin[k];
            float w = 0.0f;

            if (mel_k > left_m && mel_k < center_m) {
                w = (mel_k - left_m) / std::max(center_m - left_m, 1e-12f);
            } else if (mel_k >= center_m && mel_k < right_m) {
                w = (right_m - mel_k) / std::max(right_m - center_m, 1e-12f);
            } else {
                w = 0.0f;
            }

            if (w < 0.0f) w = 0.0f;
            fb(m, k) = w;
        }
    }
    return fb;
}

int main() {
    // ---------- Load wav (float32) ----------
    SF_INFO sfinfo{};
    SNDFILE* sndfile = sf_open("chunk0.wav", SFM_READ, &sfinfo);
    if (!sndfile) {
        std::cerr << "Error: cannot open chunk0.wav\n";
        return 1;
    }
    if (sfinfo.samplerate != SAMPLE_RATE) {
        std::cerr << "Error: expected " << SAMPLE_RATE << " Hz audio\n";
        sf_close(sndfile);
        return 1;
    }
    if (sfinfo.channels != 1) {
        std::cerr << "Error: expected mono audio (got " << sfinfo.channels << " channels)\n";
        sf_close(sndfile);
        return 1;
    }

    std::vector<float> audio(sfinfo.frames);
    sf_readf_float(sndfile, audio.data(), sfinfo.frames);
    sf_close(sndfile);

    std::cout << "Waveform shape: [1," << audio.size() << "]\n";
    std::cout.setf(std::ios::fixed); std::cout.precision(10);
    std::cout << "First 20 samples: ";
    for (int i = 0; i < 20 && i < (int)audio.size(); ++i) std::cout << audio[i] << " ";
    std::cout << "\n";

    // ---------- Global pre-emphasis (matches simple NumPy path used in HF) ----------
    if (!audio.empty()) {
        for (size_t i = audio.size() - 1; i > 0; --i) {
            audio[i] = audio[i] - PREEMPHASIS * audio[i - 1];
        }
        // audio[0] unchanged (as in many refs)
    }

    // ---------- STFT params ----------
    const int n_fft_bins = N_FFT / 2 + 1;
    const int n_frames = 1 + int((audio.size() - FRAME_LENGTH) / HOP_LENGTH);
    if (n_frames <= 0) {
        std::cerr << "Error: audio too short.\n";
        return 1;
    }

    // Hann window (periodic=False): denom = FRAME_LENGTH - 1
    std::vector<float> window(FRAME_LENGTH);
    float window_energy = 0.0f;
    for (int i = 0; i < FRAME_LENGTH; ++i) {
        float w = 0.5f - 0.5f * std::cos(2.0f * float(M_PI) * float(i) / float(FRAME_LENGTH - 1));
        window[i] = w;
        window_energy += w * w;
    }

    // STFT buffer and output
    std::vector<float> frame(N_FFT, 0.0f);
    fftwf_complex* out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n_fft_bins);
    Eigen::MatrixXf power_spec = Eigen::MatrixXf::Zero(n_fft_bins, n_frames);

    for (int f = 0; f < n_frames; ++f) {
        int start = f * HOP_LENGTH;

        // 1) copy raw frame
        std::fill(frame.begin(), frame.end(), 0.0f);
        for (int i = 0; i < FRAME_LENGTH; ++i) {
            frame[i] = audio[start + i];
        }

        // 2) remove DC offset per frame (matches HF NumPy path)
        float mean = 0.0f;
        for (int i = 0; i < FRAME_LENGTH; ++i) mean += frame[i];
        mean /= float(FRAME_LENGTH);
        for (int i = 0; i < FRAME_LENGTH; ++i) frame[i] -= mean;

        // 3) apply window
        for (int i = 0; i < FRAME_LENGTH; ++i) frame[i] *= window[i];

        // 4) FFT (float32 API)
        fftwf_plan plan = fftwf_plan_dft_r2c_1d(N_FFT, frame.data(), out, FFTW_ESTIMATE);
        fftwf_execute(plan);
        fftwf_destroy_plan(plan);

        // 5) power spectrum + normalization (window energy); single-sided correction
        for (int k = 0; k < n_fft_bins; ++k) {
            float re = out[k][0];
            float im = out[k][1];
            float p  = (re * re + im * im) / std::max(window_energy, 1e-20f);
            if (k > 0 && k < n_fft_bins - 1) p *= 2.0f;
            power_spec(k, f) = p;  // (n_fft_bins, n_frames)
        }
    }
    fftwf_free(out);

    // ---------- Mel filterbank (Kaldi scale, triangles in mel space) ----------
    Eigen::MatrixXf mel_fb = build_mel_filterbank();
    Eigen::MatrixXf mel = mel_fb * power_spec;      // (N_MELS, n_frames)

    // ---------- Natural log with floor ----------
    for (int i = 0; i < mel.rows(); ++i) {
        for (int j = 0; j < mel.cols(); ++j) {
            float v = mel(i, j);
            mel(i, j) = std::log(std::max(v, MEL_FLOOR));
        }
    }

    std::cout << "FBank raw shape: [" << mel.cols() << "," << mel.rows() << "]\n";
    std::cout << "First frame (10 bins): ";
    for (int i = 0; i < 10; ++i) std::cout << mel(i, 0) << " ";
    std::cout << "\n";

    // ---------- Pad/truncate to 1024 frames ----------
    Eigen::MatrixXf mel_fixed = Eigen::MatrixXf::Zero(N_MELS, TARGET_FRAMES);
    int copy_frames = std::min(TARGET_FRAMES, (int)mel.cols());
    mel_fixed.block(0, 0, N_MELS, copy_frames) = mel.block(0, 0, N_MELS, copy_frames);

    std::cout << "FBank padded shape: [" << mel_fixed.cols() << "," << mel_fixed.rows() << "]\n";

    // ---------- AST normalization ----------
    mel_fixed = (mel_fixed.array() - MEAN) / (STD * 2.0f);

    std::cout << "First frame after normalization (10 bins): ";
    for (int i = 0; i < 10; ++i) std::cout << mel_fixed(i, 0) << " ";
    std::cout << "\n";

    return 0;
}
