/***
  This file is part of PulseAudio.

  PulseAudio is free software; you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as published
  by the Free Software Foundation; either version 2.1 of the License,
  or (at your option) any later version.

  PulseAudio is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with PulseAudio; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
  USA.
***/

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <iomanip>

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <errno.h>

#include <pulse/simple.h>
#include <pulse/error.h>
#include <pulse/gccmacro.h>


#include "kfr/include/kfr/math.hpp"
#include "kfr/include/kfr/dft.hpp"
#include "kfr/include/kfr/dsp.hpp"
#include "kfr/include/kfr/io.hpp"

#define BUFSIZE 1024
#define FFTSIZE 64
#define NB_BINS FFTSIZE/2 // subset of the DFT bins to display/analyse -- must be at most FFTSIZE/2
//#define NB_BINS 32 // subset of the DFT bins to display/analyse

//.rate = 44100, // 8000, 22050, 11025...
//#define SAMPLERATE 44100
#define SAMPLERATE 11025

#define CALIBRATION_STEPS 10

#define DISPLAY_WIDTH 80

using namespace kfr;
using namespace std;

unsigned int nb_calibration_steps = 0;
float Amin, Amax;

float prep(float val) {

    if (nb_calibration_steps < CALIBRATION_STEPS) {
        if (val < Amin) Amin = val;
        if (val > Amax) Amax = val;
        nb_calibration_steps++;
        return 0;
    }

    return max(0, (val - Amin) / ( Amax - Amin));

}

void plothist(univector<float, FFTSIZE> data) {

    unsigned int hz_inc = SAMPLERATE / FFTSIZE;

    //size_t NBLINES = 32;
    //size_t SAMPLES_PER_LINE = FFTSIZE/NBLINES;

    cout << "\x1b[" << NB_BINS + 2 << "A\x1b[0J"; // up 10 lines, clear screen to bottom

    cout << "                 0                  25                  50                  75                 100" << endl;
    cout << "                 |-------------------|-------------------|-------------------|-------------------|" << endl;
    // The first coefficient in your array is the 0 frequency coefficient. That is basically the average power level for all frequencies. We skip it.
    // Then, we can only measure frequencies up to half the sample points. 
    unsigned int hz = 0;
    for (size_t i = 1; i < NB_BINS; i++) {
        //float avg = mean(data.slice(i * SAMPLES_PER_LINE, (i+1) * SAMPLES_PER_LINE - 1));
        cout << setfill(' ') << setw(5) << hz << " - ";
        hz += hz_inc;
        cout << setfill(' ') << setw(5) << hz << " Hz ";
        for (size_t j = 0; j < (int) (prep(data[i]) * DISPLAY_WIDTH); j++) {
            cout << "=";
        }
        cout << endl;

    }
}

        for (size_t j = 0; j < prep(data[i]); j++) {
            cout << "=";
        }
        cout << endl;

    }

}

void fft(const univector<complex<float>, BUFSIZE> samples) {


    // should we apply a window?
    //auto windowed_samples = window_hann(samples.size());

    univector<complex<float>, FFTSIZE> freq = scalar(qnan);

    dft_plan<float> dft(FFTSIZE);                      // initialize plan
    univector<u8> temp(dft.temp_size);       // allocate work buffer
    dft.execute(freq, samples, temp);        // do the actual transform
    freq = freq / FFTSIZE;
    univector<float, FFTSIZE> dB = amp_to_dB(cabs(freq));


//    println("max  = ", maxof(dB));
//    println("min  = ", minof(dB));
//    println("mean = ", mean(dB));
//    println("rms  = ", rms(dB));
//
    //println(dB);

    plothist(dB);
}

/* A simple routine calling UNIX write() in a loop */
static ssize_t loop_write(int fd, const void*data, size_t size) {
    ssize_t ret = 0;

    while (size > 0) {
        ssize_t r;

        if ((r = write(fd, data, size)) < 0)
            return r;

        if (r == 0)
            break;

        ret += r;
        data = (const uint8_t*) data + r;
        size -= (size_t) r;
    }

    return ret;
}

int main(int argc, char*argv[]) {
    /* The sample type to use */
    static const pa_sample_spec ss = {
        //.format = PA_SAMPLE_S16LE, // see https://freedesktop.org/software/pulseaudio/doxygen/sample_8h.html#a3c622fc51f4fc6ebfdcc7b454ac9c05f
        .format = PA_SAMPLE_FLOAT32LE, // see https://freedesktop.org/software/pulseaudio/doxygen/sample_8h.html#a3c622fc51f4fc6ebfdcc7b454ac9c05f
        .rate = SAMPLERATE,
        .channels = 1
    };
    pa_simple *s = NULL;
    int ret = 1;
    int error;

    /* Create the recording stream */
    if (!(s = pa_simple_new(NULL, argv[0], PA_STREAM_RECORD, NULL, "record", &ss, NULL, NULL, &error))) {
        fprintf(stderr, __FILE__": pa_simple_new() failed: %s\n", pa_strerror(error));
        goto finish;
    }

    for (;;) {
        float buf[BUFSIZE];

        /* Record some data ... */
        if (pa_simple_read(s, buf, sizeof(buf), &error) < 0) {
            fprintf(stderr, __FILE__": pa_simple_read() failed: %s\n", pa_strerror(error));
            goto finish;
        }

        /* And write it to STDOUT */
        //if (loop_write(STDOUT_FILENO, buf, sizeof(buf)) != sizeof(buf)) {
        //    fprintf(stderr, __FILE__": write() failed: %s\n", strerror(errno));
        //    goto finish;
        //}
        univector<complex<float>, BUFSIZE> samples;
        for(size_t i = 0; i < BUFSIZE; i++) {
            samples[i] = make_complex(buf[i]);
        }
        fft(samples);
    }

    ret = 0;

finish:

    if (s)
        pa_simple_free(s);

    return ret;
}
