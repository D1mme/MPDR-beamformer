# Minimum Power Distortionless Response Beamformer
This code implements a near-field minimum power distortionless response (MPDR) beamformer.

As input, the beamformer takes the target location and the microphone locations. The power spectral density matrices corresponding to the received signals are estimated. 
Using this matrix, the MPDR beamformer tries to preserve the signal from the target location while discarding the noise signal. 
If a voice activity detector is used, the MPDR can be changed to a minimum variance distortionless response (MVDR) beamformer.

Three examples are provided. 
- Example 1: An anechoic scenario with a single noise source
- Example 2: A reverberant scenario with a single noise source
- Example 3: A reverberant scenario with multiple noise sources

The examples make use of the [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) from E. Habets (MIT license). You might need to compile this for your system.
The sound excerpt is taken from the movie ['Sprite Fight'](https://studio.blender.org/films/sprite-fright/) by Blender Studio (Creative Commons Attribution 1.0 License). 


 
