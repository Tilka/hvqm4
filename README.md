# HVQM4 1.3/1.5 (*.h4m) video decoder

This code is in a bit of a state right now and needs someone to tidy it up and integrate it into FFmpeg.
I probably won't get around to it.
You can tell it's a hurried holiday project because the file is still called h4m_audio_decode.c.
However, it is bit-accurate for both 1.3 and 1.5, and ABI compatible with 1.5.
The ABI compatibility can probably be dropped now that it works.
Feel free to ask questions as long as it's not "how do I run this on Windows".
The GameCube game "PK: Out of the Shadows" appears to be the only game that links the debug version of the HVQM4 SDK.
This is useful because it doesn't inline any functions.

Ignore the game "Frogger Beyond", it uses some strange intermediate version of the SDK.
Its videos should decode just fine though.

Also ignore the audio part, I believe vgmstream has some fixes that this version of h4m_audio_decode doesn't have.

## Testing
yoloader is a basic PowerPC ELF loader I wrote for testing against the original Hudson Soft decoder.
I used [buildroot](https://buildroot.org/) to build a crosscompiler toolchain for 32-bit PowerPC Linux targets and then ran the decoder in QEMU usermode emulation.
The HVQM4 SDK doesn't require any MMIO and doesn't do anything weird ABI-wise, so it just works.
You can load game executables that have function symbols and relocation information (such as the PK executable) and replace individual functions.
If function foo() calls bar(), you can test foo() in isolation by changing foo() to call pbar() instead.
Decoding of videos in portrait orientation is untested (there are no games that use it).

## Credits
This code is based on the HVQM4 1.3/1.5 audio decoder by hcs.

## License
LGPLv2 or later
