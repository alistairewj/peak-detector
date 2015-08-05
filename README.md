# peak-detector
This repository contains code for multimodal R-peak detection code. The R-peak is the prominent portion of the QRS complex - a regularly occuring pattern on an electrocardiogram (ECG) that corresponds to a heart beat. The software here aims to precisely locate the R-peak using not only the ECG, but also the arterial blood pressure (ABP) waveform, the photoplethysmogram (PPG) and/or the stroke volume (SV).

Briefly, the technique aims to fuse the signals based off an estimate of signal quality known as a signal quality index (SQI). An SQI is estimated for each type of signal (e.g. ECG, ABP, etc) and the peak detections for each signal are fused if and only if their SQI is above a threshold. Additional code is implemented which accounts for delay present on signals not directly measuring the heart beat (e.g. the pulsatile waveform on the PPG usually occurs much later than the corresponding QRS complex in the ECG).

# Acknowledgement

The code was developed for the PhysioNet/Computing in Cardiology 2014 challenge on multimodal peak detection. A conference publication [1] describing the principles of the technique is available here: http://cinc.org/archives/2014/pdf/0281.pdf
The current implementation, as of June 1st, 2015 is described in a [Physiological Measurement special issue article](http://stacks.iop.org/0967-3334/36/1665) [2].


[1] Johnson, A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2014).  R-peak estimation using multimodal  lead switching, Computing in Cardiology Conference (CinC), 2014, Vol. 41, pp. 281-284.

[2] Johnson,  A. E. W., Behar, J., Andreotti, F., Clifford, G. D. and Oster, J. (2015). Multimodal heart beat detection using signal quality indices, Physiological Measurement 36 (2015): 1665-1677.


