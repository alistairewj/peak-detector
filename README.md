# peak-detector
This repository contains MATLAB code for multimodal R-peak detection code. The R-peak is the prominent portion of the QRS complex - a regularly occuring pattern on an electrocardiogram (ECG) that corresponds to a heart beat. The software here aims to precisely locate the R-peak using not only the ECG, but also the arterial blood pressure (ABP) waveform, the photoplethysmogram (PPG) and/or the stroke volume (SV).

Briefly, the technique aims to fuse the signals based off an estimate of signal quality known as a signal quality index (SQI). An SQI is estimated for each type of signal (e.g. ECG, ABP, etc) and the peak detections for each signal are fused if and only if their SQI is above a threshold. Additional code is implemented which accounts for delay present on signals not directly measuring the heart beat (e.g. the pulsatile waveform on the PPG usually occurs much later than the corresponding QRS complex in the ECG).

# Requirements

All code within this repository requires the WFDB toolbox:
https://github.com/ikarosilva/wfdb-app-toolbox

After downloading the above package, ensure the subfolder 'mcode' is on your path in order to run these functions.

# Using this repository

There are four main functions:

* **detect.m** - The primary function, which runs the algorithm on a WFDB-readable record
* **detect_matlab.m** - Same as the above, but runs the algorithm on a MATLAB data matrix with associated header and sampling frequency
* **detect_sqi.m** - A subroutine of detect.m - this runs the SQI switching approach
* **detect_regularity.m** - A subroutine of detect.m - this runs the regularity switching approach

There are two examples provided in the examples subfolder to help you get started.

# More information about functions

## detect_matlab

The function `detect_matlab` has been provided to perform peak detection directly on data stored in MATLAB. The function can be called as follows:

```matlab
[ qrs, sqi ] = detect_matlab(DATA, HEADER, FS, OPT)
```

This function implements the SQI switching method which had the highest performance in Phase 3 of the Physionet/Computing in Cardiology 2014 challenge on multimodal peak detection. The inputs are:

* DATA - An NxM matrix of N samples with M signals. Data should be input in physical units, otherwise the estimates of signal quality will fail.
* HEADER - A 1xM cell array of strings which provides signal information for each signal in DATA. Please use recognizable signal names: 'ecg' will be recognized as an electrocardiogram, but 'signal 1' will not.
* FS - A 1x1 scalar which represents the sampling frequency at which the data is stored.
* OPT_INPUT - This is a structure used to define a variety of parameters for each peak detector (e.g. the ECG peak detector, the ABP peak detector, etc).

More information on the options structure is provided later in this readme.

## detect.m

The main function for the peak detection algorithm is [detect.m](https://github.com/alistairewj/peak-detector/blob/master/detect.m). The syntax is as follows:

```matlab
[ qrs_final, sqi_ecg, ann_jqrs, ann_gqrs ] = detect(recordName, alg, opt)
```

The four inputs are:

* `recordName` - this is the name of the file to be loaded in with the WFDB toolbox. If you are unsure what the WFDB toolbox is, please see [here](http://physionet.org/physiotools/matlab/wfdb-app-matlab/). There should be two files: one `.dat` and one `.hea`.
* ALG - This string defines the type of fusion to do across leads. `'sqi'` specifies that signal quality should be used to fuse the leads. `'regularity'` specifies that regularity should be used. The default (for any text, including `'sqi'`) is to use `'sqi'`.
* opt - This is a structure used to define a variety of parameters for each peak detector (e.g. the ECG peak detector, the ABP peak detector, etc).

More information on the options structure is provided below.

## Options structure

opt_input has the following fields (all set to sensible defaults if not otherwise specified):

    SIZE_WIND - integer (seconds) - the size of the window to use for evaluating signal quality
    LG_MED - integer (measurements) - take the median SQI using X nearby values, used only for the ECG SQI, so if LG_MED = 3, we take the median of the 3 prior and 3 posterior windows
    REG_WIN - integer (seconds) - how frequently to check the SQI for switching, e.g. if REG_WIN=1 then we can switch the signal used every second
    THR - double (seconds) - the maximum window allowed to consider to peaks "matched" in the calculation of the ECG SQI
    SQI_THR - double ([0-1]) - minimum SQI required to use this signal
    USE_PACING - binary (0 or 1) - whether to crudely detect and compensate for pacing signals
    ABPMethod - string - ABP peak detection method, either 'wabp' or 'delineator'
    SIMPLEMODE - binary (0 or 1) - simple mode uses only the first occuring ABP and first occuring ECG signal
    JQRS_THRESH - double - threshold to determine a peak in the ECG peak detector 'jqrs'
    JQRS_REFRAC - double (seconds) - refractory period used in the ECG peak detector 'jqrs'
    JQRS_INTWIN_SZ - integer
    JQRS_WINDOW - integer - window size used for repeated applications of 'jqrs'
    DELAY - This is the method to determine the delay between the ECG R-peak and any associated pulse in another waveform. The default, `'map'`, finds a 1 minute segment with alternating ECG R-peaks and pulses, and calculates the average delay for this segment. The alternative method, `'crosscorr'`, uses cross correlation to determine the delay. The single delay value for both methods is used across the *entire* signal. This was valid for the 10 minute signals that the code was developed on: it may not be valid on longer records.

# Acknowledgement

The algorithm was described in [Physiological Measurement special issue article](http://stacks.iop.org/0967-3334/36/1665) [1], and was an extension to an earlier conference publication [2].

If you use this code in your research, we would be grateful if you acknowledged it with a citation to [1].

[1] Johnson,  A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2015). Multimodal heart beat detection using signal quality indices, Physiological Measurement 36 (2015): 1665-1677.

[2] Johnson, A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2014).  R-peak estimation using multimodal  lead switching, Computing in Cardiology Conference (CinC), 2014, Vol. 41, pp. 281-284.
