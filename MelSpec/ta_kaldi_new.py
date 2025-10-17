


window_function = _feature_window_function(window_type = "hanning", window_size = 400, blackman_coeff = 0.42, device, dtype).unsqueeze(
        0
    )  # size (1, window_size)

   
   




### This is the function called from the feature extractor.py

def fbank(
    waveform: Tensor,
    blackman_coeff: float = 0.42,
    channel: int = -1,
    dither: float = 0.0,
    energy_floor: float = 1.0,
    frame_length: float = 25.0,                     # 25 * fs / 1000.  ( 400 )
    frame_shift: float = 10.0,                      # 10 * fs / 1000   ( 160 )
    high_freq: float = 0.0,
    htk_compat: bool = False,
    low_freq: float = 20.0,
    min_duration: float = 0.0,
    num_mel_bins: int = 23,                         # 128
    preemphasis_coefficient: float = 0.97,
    raw_energy: bool = True,
    remove_dc_offset: bool = True,
    round_to_power_of_two: bool = True,
    sample_frequency: float = 16000.0,              # 16000.0
    snip_edges: bool = True,
    subtract_mean: bool = False,
    use_energy: bool = False,
    use_log_fbank: bool = True,                     # log 
    use_power: bool = True,                         # true
    vtln_high: float = -500.0,
    vtln_low: float = 100.0,
    vtln_warp: float = 1.0,
    window_type: str = POVEY,                        # "hanning"
) -> Tensor:
    r"""Create a fbank from a raw audio signal. This matches the input/output of Kaldi's
    compute-fbank-feats.

    Args:
        waveform (Tensor): Tensor of audio of size (c, n) where c is in the range [0,2)
        blackman_coeff (float, optional): Constant coefficient for generalized Blackman window. (Default: ``0.42``)
        channel (int, optional): Channel to extract (-1 -> expect mono, 0 -> left, 1 -> right) (Default: ``-1``)
        dither (float, optional): Dithering constant (0.0 means no dither). If you turn this off, you should set
            the energy_floor option, e.g. to 1.0 or 0.1 (Default: ``0.0``)
        energy_floor (float, optional): Floor on energy (absolute, not relative) in Spectrogram computation.  Caution:
            this floor is applied to the zeroth component, representing the total signal energy.  The floor on the
            individual spectrogram elements is fixed at std::numeric_limits<float>::epsilon(). (Default: ``1.0``)
        frame_length (float, optional): Frame length in milliseconds (Default: ``25.0``)
        frame_shift (float, optional): Frame shift in milliseconds (Default: ``10.0``)
        high_freq (float, optional): High cutoff frequency for mel bins (if <= 0, offset from Nyquist)
         (Default: ``0.0``)
        htk_compat (bool, optional): If true, put energy last.  Warning: not sufficient to get HTK compatible features
         (need to change other parameters). (Default: ``False``)
        low_freq (float, optional): Low cutoff frequency for mel bins (Default: ``20.0``)
        min_duration (float, optional): Minimum duration of segments to process (in seconds). (Default: ``0.0``)
        num_mel_bins (int, optional): Number of triangular mel-frequency bins (Default: ``23``)
        preemphasis_coefficient (float, optional): Coefficient for use in signal preemphasis (Default: ``0.97``)
        raw_energy (bool, optional): If True, compute energy before preemphasis and windowing (Default: ``True``)
        remove_dc_offset (bool, optional): Subtract mean from waveform on each frame (Default: ``True``)
        round_to_power_of_two (bool, optional): If True, round window size to power of two by zero-padding input
            to FFT. (Default: ``True``)
        sample_frequency (float, optional): Waveform data sample frequency (must match the waveform file, if
            specified there) (Default: ``16000.0``)
        snip_edges (bool, optional): If True, end effects will be handled by outputting only frames that completely fit
            in the file, and the number of frames depends on the frame_length.  If False, the number of frames
            depends only on the frame_shift, and we reflect the data at the ends. (Default: ``True``)
        subtract_mean (bool, optional): Subtract mean of each feature file [CMS]; not recommended to do
            it this way.  (Default: ``False``)
        use_energy (bool, optional): Add an extra dimension with energy to the FBANK output. (Default: ``False``)
        use_log_fbank (bool, optional):If true, produce log-filterbank, else produce linear. (Default: ``True``)
        use_power (bool, optional): If true, use power, else use magnitude. (Default: ``True``)
        vtln_high (float, optional): High inflection point in piecewise linear VTLN warping function (if
            negative, offset from high-mel-freq (Default: ``-500.0``)
        vtln_low (float, optional): Low inflection point in piecewise linear VTLN warping function (Default: ``100.0``)
        vtln_warp (float, optional): Vtln warp factor (only applicable if vtln_map not specified) (Default: ``1.0``)
        window_type (str, optional): Type of window ('hamming'|'hanning'|'povey'|'rectangular'|'blackman')
         (Default: ``'povey'``)

    Returns:
        Tensor: A fbank identical to what Kaldi would output. The shape is (m, ``num_mel_bins + use_energy``)
        where m is calculated in _get_strided
    """
   
    #device, dtype = waveform.device, waveform.dtype
    # dtypt = float32

    #waveform, window_shift, window_size, padded_window_size = _get_waveform_and_window_properties(
    #    waveform, channel, sample_frequency, frame_shift, frame_length, round_to_power_of_two, preemphasis_coefficient
    #)

    """
    # waveform in mono, 
    # hoplen = 160
    # winsize = 400 
    # nfft = 512
    """
    

    #if len(waveform) < min_duration * sample_frequency:
    #    # signal is too short
    #    return torch.empty(0, device=device, dtype=dtype)

    # strided_input, size (m, padded_window_size) and signal_log_energy, size (m)
    
    strided_input, signal_log_energy = _get_window(
        waveform,
        padded_window_size,
        window_size,
        window_shift,
        window_type,
        blackman_coeff,
        snip_edges,
        raw_energy,
        energy_floor,
        dither,
        remove_dc_offset,
        preemphasis_coefficient,
    )

    # size (m, padded_window_size // 2 + 1)
    spectrum = torch.fft.rfft(strided_input).abs()
    if use_power:
        spectrum = spectrum.pow(2.0)


    # gets the mel filterbank
    # pads them at right - is this DC?
    # multiplies by the spectrum
        # sets zeros to epsilon    

    # size (num_mel_bins, padded_window_size // 2)
    mel_energies, _ = get_mel_banks(
        num_mel_bins, padded_window_size, sample_frequency, low_freq, high_freq, vtln_low, vtln_high, vtln_warp
    )
    mel_energies = mel_energies.to(device=device, dtype=dtype)


    # pad right column with zeros and add dimension, size (num_mel_bins, padded_window_size // 2 + 1)
    mel_energies = torch.nn.functional.pad(mel_energies, (0, 1), mode="constant", value=0)


    # sum with mel fiterbanks over the power spectrum, size (m, num_mel_bins)
    mel_energies = torch.mm(spectrum, mel_energies.T)
    
    if use_log_fbank:
        # avoid log of zero (which should be prevented anyway by dithering)
        mel_energies = torch.max(mel_energies, _get_epsilon(device, dtype)).log()

    # if use_energy then add it as the last column for htk_compat == true else first column
    #if use_energy:
    #    signal_log_energy = signal_log_energy.unsqueeze(1)  # size (m, 1)
    #    # returns size (m, num_mel_bins + 1)
    #    if htk_compat:
    #        mel_energies = torch.cat((mel_energies, signal_log_energy), dim=1)
    #    else:
    #        mel_energies = torch.cat((signal_log_energy, mel_energies), dim=1)


    #mel_energies = _subtract_column_mean(mel_energies, subtract_mean)
    #return mel_energies





def _get_window(
    waveform: Tensor,
    padded_window_size: int,                 # 512
    window_size: int,                        # 400
    window_shift: int,                       # 160
    window_type: str,                        # hanning
    blackman_coeff: float,                   # 0.42
    snip_edges: bool,                        # True
    raw_energy: bool,                        # True
    energy_floor: float,                     # 1.0
    dither: float,                           # 0.0
    remove_dc_offset: bool,                  # True
    preemphasis_coefficient: float,          #0.97
) -> Tuple[Tensor, Tensor]:
    r"""Gets a window and its log energy

    Returns:
        (Tensor, Tensor): strided_input of size (m, ``padded_window_size``) and signal_log_energy of size (m)
    """
    device, dtype = waveform.device, waveform.dtype
    epsilon = _get_epsilon(device, dtype)

    # size (m, window_size)
    strided_input = _get_strided(waveform, window_size, window_shift, snip_edges)

    #if dither != 0.0:
    #    rand_gauss = torch.randn(strided_input.shape, device=device, dtype=dtype)
    #    strided_input = strided_input + rand_gauss * dither

    if remove_dc_offset:
        # Subtract each row/frame by its mean
        row_means = torch.mean(strided_input, dim=1).unsqueeze(1)  # size (m, 1)
        strided_input = strided_input - row_means

    if raw_energy:
        # Compute the log energy of each row/frame before applying preemphasis and
        # window function
        signal_log_energy = _get_log_energy(strided_input, epsilon, energy_floor)  # size (m)

    if preemphasis_coefficient != 0.0:
        # strided_input[i,j] -= preemphasis_coefficient * strided_input[i, max(0, j-1)] for all i,j
        offset_strided_input = torch.nn.functional.pad(strided_input.unsqueeze(0), (1, 0), mode="replicate").squeeze(
            0
        )  # size (m, window_size + 1)
        strided_input = strided_input - preemphasis_coefficient * offset_strided_input[:, :-1]


    window_function = torch.hann_window(400, periodic=False, device=device, dtype=dtype)
    strided_input = strided_input * window_function  # size (m, window_size)


    if padded_window_size != window_size:
        # Pad columns with zero until we reach size (m, padded_window_size)
        padding_right = padded_window_size - window_size
        strided_input = torch.nn.functional.pad(
            strided_input.unsqueeze(0), (0, padding_right), mode="constant", value=0
        ).squeeze(0)

    # Compute energy after window function (not the raw one)
    #if not raw_energy:
    #    signal_log_energy = _get_log_energy(strided_input, epsilon, energy_floor)  # size (m)

    return strided_input, signal_log_energy
