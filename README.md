# Minimum Power Distortionless Response Beamformer
This code is the matlab implementation of  a near-field minimum power distortionless response (MPDR) beamformer. 
Note that, in literature, the MPDR beamformer is often incorrectly referred to as the MVDR beamformer.

As input, the beamformer takes the target location and the microphone locations. The power spectral density matrix corresponding to the received signal is estimated. 
Using this matrix, the MPDR beamformer tries to preserve the signal from the target location while discarding signals from other directions. 
If a voice activity detector is used, the MPDR can be changed to a minimum variance distortionless response (MVDR) beamformer.

Four examples are provided. 
- Example 1: An anechoic scenario with a single noise source, MPDR beamforming
- Example 2: A reverberant scenario with a single noise source, MPDR beamforming
- Example 3: A reverberant scenario with multiple noise sources, MPDR beamforming
- Example 4: An anechoic scenario with a single noise source, MVDR beamforming with perfect voice activity detection

For each of the examples, three output audio files can be found in `Example_Audio`: 
- The clean (target) signal. Note that it includes a small delay before starting
- The noisy signal. The noise level is not necessarily the same throughout the different examples. This signal is the signal measured by the reference microphone. The reference microphone is the microphone nearest to the target.
- The enhanced signal. I.e. the output of the microphone beamformer. 

## Licensing
The examples make use of the [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) from E. Habets (MIT license). You might need to compile this for your system.
The sound excerpt is taken from the movie ['Sprite Fight'](https://studio.blender.org/films/sprite-fright/) by Blender Studio (Creative Commons Attribution 1.0 License). 

## Citation information:
We used a variant of this implementation in our ICASSP paper "Loudspeaker Beamforming to Enhance Speech Recognition Performance of Voice Driven Applications". If you use this implementation please consider citing us. Citation information:
```
@INPROCEEDINGS{10889702,
  author={de Groot, D. and Karslioglu, B. and Scharenborg, O. and Martinez, J.},
  booktitle={ICASSP 2025 - 2025 IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP)},
  title={Loudspeaker Beamforming to Enhance Speech Recognition Performance of Voice Driven Applications},
  year={2025},
  volume={},
  number={},
  pages={1-5},
  keywords={Loudspeakers;Performance evaluation;Acoustic distortion;Array signal processing;Signal processing algorithms;Acoustics;Robustness;Distortion measurement;Speech processing;Automatic speech recognition;Spotforming;beamforming;speech recognition}
  doi={10.1109/ICASSP49660.2025.10889702}}
```

