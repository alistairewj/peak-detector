# peak-detector
This repository contains MATLAB code for multimodal R-peak detection code. The R-peak is the prominent portion of the QRS complex - a regularly occuring pattern on an electrocardiogram (ECG) that corresponds to a heart beat. The software here aims to precisely locate the R-peak using not only the ECG, but also the arterial blood pressure (ABP) waveform, the photoplethysmogram (PPG) and/or the stroke volume (SV).

Briefly, the technique aims to fuse the signals based off an estimate of signal quality known as a signal quality index (SQI). An SQI is estimated for each type of signal (e.g. ECG, ABP, etc) and the peak detections for each signal are fused if and only if their SQI is above a threshold. Additional code is implemented which accounts for delay present on signals not directly measuring the heart beat (e.g. the pulsatile waveform on the PPG usually occurs much later than the corresponding QRS complex in the ECG).

# How to use this repository

The main function for the peak detection algorithm is [detect.m](https://github.com/alistairewj/peak-detector/blob/master/detect.m). The syntax is as follows:

```matlab
[ qrs_final, sqi_ecg, ann_jqrs, ann_gqrs ] = detect(recordName, FUSEALG, DELAYALG, opt_input)
```

The four inputs are:

* `recordName` - this is the name of the file to be loaded in with the WFDB toolbox. If you are unsure what the WFDB toolbox is, please see [here](http://physionet.org/physiotools/matlab/wfdb-app-matlab/). There should be two files: one `.dat` and one `.hea`.
* FUSEALG - This string defines the type of fusion to do across leads. `'sqi'` specifies that signal quality should be used to fuse the leads. `'regularity'` specifies that regularity should be used. The default (for any text, including `'sqi'`) is to use `'sqi'`.
* DELAY - This is the method to determine the delay between the ECG R-peak and any associated pulse in another waveform. The default, `'map'`, finds a 1 minute segment with alternating ECG R-peaks and pulses, and calculates the average delay for this segment. This delay is used across the *entire* signal. This was valid for the 10 minute signals that the code was developed on: it may not be valid on longer records. The alternative method, `'crosscorr'`, uses cross correlation to determine the delay.
* opt_input - This is a structure used to define a variety of parameters for each peak detector (e.g. the ECG peak detector, the ABP peak detector, etc).

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

# Acknowledgement

The code was developed for the PhysioNet/Computing in Cardiology 2014 challenge on multimodal peak detection. A conference publication [1] describing the principles of the technique is available here: http://cinc.org/archives/2014/pdf/0281.pdf

The current implementation, as of June 1st, 2015 is described in a [Physiological Measurement special issue article](http://stacks.iop.org/0967-3334/36/1665) [2].

If you use this code in your research, we would be grateful if you acknowledged it with a citation to [2].

[1] Johnson, A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2014).  R-peak estimation using multimodal  lead switching, Computing in Cardiology Conference (CinC), 2014, Vol. 41, pp. 281-284.

[2] Johnson,  A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2015). Multimodal heart beat detection using signal quality indices, Physiological Measurement 36 (2015): 1665-1677.


