import os
import soundfile as sf
import torchaudio
import torch
import numpy as np
from transformers import AutoFeatureExtractor

# -------------------------------
# 1) Paths + Settings
# -------------------------------
CACHE_DIR = "cache"
EXTRACTOR_PATH = os.path.join(CACHE_DIR, "extractor")  # local extractor cache

AUDIO_FOLDER = "unseen_files"
OUT_DIR = "features"
VALID_EXTS = (".wav", ".mp3", ".flac", ".aiff", ".ogg", ".m4a")
SAMPLE_RATE = 16000

os.makedirs(OUT_DIR, exist_ok=True)

# -------------------------------
# 2) Load Feature Extractor from cache
# -------------------------------
feature_extractor = AutoFeatureExtractor.from_pretrained(EXTRACTOR_PATH)
print("📂 Feature extractor loaded from cache")

# -------------------------------
# 3) Audio Loader
# -------------------------------
def load_audio(path, target_sr=16000):     # problem with mono and downsampling
    wav, sr = sf.read(path)
    if wav.ndim > 1:  # stereo → mono
        wav = np.mean(wav, axis=1)
    if sr != target_sr:
        wav = torchaudio.functional.resample(
            torch.tensor(wav), sr, target_sr
        ).numpy()
        sr = target_sr
    return wav, sr

# -------------------------------
# 4) Extract + Save Features
# -------------------------------
def extract_and_save(file_path):
    try:
        wav, sr = load_audio(file_path, SAMPLE_RATE)

        # Extract features (AST input log-mel patches)
        features = feature_extractor(
            wav,
            sampling_rate=sr,
            return_tensors="np"   # return NumPy arrays
        )["input_values"]


        # Save to .npy file
        fname = os.path.splitext(os.path.basename(file_path))[0] + "_features.npy"
        out_path = os.path.join(OUT_DIR, fname)
        np.save(out_path, features)

        print(f"✅ Saved features: {out_path}")
    except Exception as e:
        print(f"❌ Error with {file_path}: {e}")


def extract_features(file_path):

    wav, sr = load_audio(file_path, SAMPLE_RATE)

    # Extract features (AST input log-mel patches)
    features = feature_extractor(
        wav,
        sampling_rate=sr,
        return_tensors="np"   # return NumPy arrays
    )["input_values"]

    return features
    

"""
# -------------------------------
# 5) Run on all files
# -------------------------------
files = [f for f in os.listdir(AUDIO_FOLDER) if f.lower().endswith(VALID_EXTS)]
print(f"🎶 Found {len(files)} audio files")

for f in files:
    file_path = os.path.join(AUDIO_FOLDER, f)
    extract_and_save(file_path)
"""


# OUR VERSION :
    # load audio with librosa ( using the resample there and Mono )
    # save that to a .wav file (scipy ? )

    # load the wav in python librosa and do the (mel) spectrogram using the feature_extractor
    
    # call the CLI objC function
    # load the ObjC outputs from file

    # Compare these outputs with the ObjC


# checks - will need to save the windowed (processed) data

# checks - will need to check the window functions - hanning 400

# checks - will need to check the outputs of the fft

# checks - will need to check the mel filterbanks