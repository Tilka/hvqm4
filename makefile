native:
	clang -Wall -Wextra -g -fsanitize=address -fno-omit-frame-pointer -m32 h4m_audio_decode.c -o h4m_audio_decode -DNATIVE=1
	rm -f *.ppm
	./h4m_audio_decode samples/LOGOS.h4m foo.wav

emu:
	toolchain/bin/powerpc-linux-gcc -Wall -Wextra -g -fno-omit-frame-pointer -static h4m_audio_decode.c -o h4m_audio_decode
	rm -f *.ppm
	qemu-ppc h4m_audio_decode samples/LOGOS.h4m foo.wav

clean:
	rm -f h4m_audio_decode *.ppm

.PHONY: all clean
