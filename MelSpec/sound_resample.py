import os
import librosa
import soundfile as sf
import numpy as np

# ==== CONFIG ====
INPUT_ROOT = "unseensample"         # Root folder with your subfolders of audio
OUTPUT_ROOT = "unseen_files"   # Where processed audio will be saved
TARGET_SR = 16000           # 16kHz (NSynth standard)
TARGET_DURATION = 10        # seconds
TARGET_SAMPLES = TARGET_SR * TARGET_DURATION
SUPPORTED_FORMATS = (".wav", ".mp3", ".aiff")

# Ensure output root exists
os.makedirs(OUTPUT_ROOT, exist_ok=True)


def process_audio(file_path, out_dir):
    """Load, process, and save audio to exactly 10 seconds at 16kHz."""
    try:
        # Load audio file (mono, resampled)
        y, sr = librosa.load(file_path, sr=TARGET_SR, mono=True)
    except Exception as e:
        print(f"❌ Error loading {file_path}: {e}")
        return

    # Skip empty files
    if y is None or len(y) == 0:
        print(f"⚠️ Skipping empty or unreadable file: {file_path}")
        return

    base_name = os.path.splitext(os.path.basename(file_path))[0]

    # If shorter than 10 sec → repeat until 10 sec
    if len(y) < TARGET_SAMPLES:
        repeats = int(np.ceil(TARGET_SAMPLES / max(len(y), 1)))
        padded = np.tile(y, repeats)[:TARGET_SAMPLES]
        out_name = f"{base_name}_short.wav"
        sf.write(os.path.join(out_dir, out_name), padded, TARGET_SR)
        print(f"✅ Processed short file → {out_name}")

    else:
        # Split into 10-sec chunks
        num_chunks = len(y) // TARGET_SAMPLES
        for i in range(num_chunks):
            chunk = y[i*TARGET_SAMPLES:(i+1)*TARGET_SAMPLES]
            out_name = f"{base_name}_chunk{i}.wav"
            sf.write(os.path.join(out_dir, out_name), chunk, TARGET_SR)

        # Handle remainder (if not exactly divisible)
        remainder = y[num_chunks*TARGET_SAMPLES:]
        if len(remainder) > 0:
            repeats = int(np.ceil(TARGET_SAMPLES / len(remainder)))
            padded = np.tile(remainder, repeats)[:TARGET_SAMPLES]
            out_name = f"{base_name}_chunk{num_chunks}_rem.wav"
            sf.write(os.path.join(out_dir, out_name), padded, TARGET_SR)
        print(f"✅ Processed long file → {base_name} into {num_chunks+1} chunks")


def main():
    """Walk through all subfolders and process audio files."""
    for root, _, files in os.walk(INPUT_ROOT):
        for file in files:
            if file.lower().endswith(SUPPORTED_FORMATS):
                file_path = os.path.join(root, file)

                # Mirror subfolder structure in output
                rel_dir = os.path.relpath(root, INPUT_ROOT)
                out_dir = os.path.join(OUTPUT_ROOT, rel_dir)
                os.makedirs(out_dir, exist_ok=True)

                process_audio(file_path, out_dir)


if __name__ == "__main__":
    main()
