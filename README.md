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

The examples make use of the [room-impulse response generator](https://www.audiolabs-erlangen.de/fau/professor/habets/software/rir-generator) from E. Habets (MIT license). You might need to compile this for your system.
The sound excerpt is taken from the movie ['Sprite Fight'](https://studio.blender.org/films/sprite-fright/) by Blender Studio (Creative Commons Attribution 1.0 License). 


 
