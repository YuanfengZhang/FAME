// ============================================================================
// zstdstream, C++ iostream classes wrapping the zstd compression library.
// Based on gzstream by Deepak Bandyopadhyay, Lutz Kettner
// Adapted for zstd by FAME project
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// ============================================================================
//
// File          : zstdstream.C
// Author(s)     : Based on gzstream, adapted for zstd
//
// Standard streambuf implementation following Nicolai Josuttis, "The
// Standard C++ Library".
// ============================================================================

#include "zstdstream.h"
#include <iostream>
#include <cstring>  // for memcpy

#ifdef ZSTDSTREAM_NAMESPACE
namespace ZSTDSTREAM_NAMESPACE {
#endif

// ----------------------------------------------------------------------------
// Internal classes to implement zstdstream. See header file for user classes.
// ----------------------------------------------------------------------------

// --------------------------------------
// class zstdstreambuf:
// --------------------------------------

zstdstreambuf* zstdstreambuf::open( const char* name, int open_mode) {
    if ( is_open())
        return (zstdstreambuf*)0;
    mode = open_mode;
    // only support input mode for decompression
    if ((mode & std::ios::ate) || (mode & std::ios::app)
        || (mode & std::ios::out))
        return (zstdstreambuf*)0;

    // Open raw file
    file = fopen( name, "rb");
    if (file == nullptr)
        return (zstdstreambuf*)0;

    // Create zstd decompression stream
    dstream = ZSTD_createDStream();
    if (dstream == nullptr) {
        fclose(file);
        file = nullptr;
        return (zstdstreambuf*)0;
    }

    // Initialize decompression stream
    size_t initResult = ZSTD_initDStream(dstream);
    if (ZSTD_isError(initResult)) {
        ZSTD_freeDStream(dstream);
        dstream = nullptr;
        fclose(file);
        file = nullptr;
        return (zstdstreambuf*)0;
    }

    // Allocate input buffer (use recommended size)
    size_t inBufSize = ZSTD_DStreamInSize();
    inBuf = new char[inBufSize];

    // Initialize buffers
    input.src = inBuf;
    input.size = 0;
    input.pos = 0;

    output.dst = buffer;
    output.size = bufferSize;
    output.pos = 0;

    opened = 1;
    return this;
}

zstdstreambuf* zstdstreambuf::close() {
    if ( is_open()) {
        sync();
        opened = 0;

        // Free zstd resources
        if (dstream != nullptr) {
            ZSTD_freeDStream(dstream);
            dstream = nullptr;
        }

        // Free input buffer
        if (inBuf != nullptr) {
            delete[] inBuf;
            inBuf = nullptr;
        }

        // Close file
        if (file != nullptr) {
            fclose(file);
            file = nullptr;
        }

        return this;
    }
    return (zstdstreambuf*)0;
}

int zstdstreambuf::underflow() { // used for input buffer only
    if ( gptr() && ( gptr() < egptr()))
        return * reinterpret_cast<unsigned char *>( gptr());

    if ( ! (mode & std::ios::in) || ! opened)
        return EOF;

    // Josuttis' implementation of inbuf
    int n_putback = gptr() - eback();
    if ( n_putback > 4)
        n_putback = 4;
    memcpy( buffer + (4 - n_putback), gptr() - n_putback, n_putback);

    // Try to decompress more data
    size_t totalDecompressed = 0;
    bool needMoreInput = true;
    bool reachedEOF = false;

    // Loop until we have some decompressed data or hit EOF
    while (totalDecompressed == 0 && !reachedEOF) {
        // If we've consumed all input, read more from file
        if (input.pos >= input.size) {
            size_t inBufSize = ZSTD_DStreamInSize();
            size_t bytesRead = fread(inBuf, 1, inBufSize, file);

            if (bytesRead == 0) {
                // EOF reached - check if there's remaining data in decompressor
                if (input.pos < input.size) {
                    // Still have data to process, continue
                } else {
                    reachedEOF = true;
                    break;
                }
            } else {
                input.size = bytesRead;
                input.pos = 0;
            }
        }

        // Set output buffer position
        output.dst = buffer + 4 + totalDecompressed;
        output.size = bufferSize - 4 - totalDecompressed;
        output.pos = 0;

        // Decompress a chunk
        size_t result = ZSTD_decompressStream(dstream, &output, &input);

        if (ZSTD_isError(result)) {
            // Decompression error
            return EOF;
        }

        totalDecompressed += output.pos;

        // If result == 0, entire frame was consumed
        // May need to read next frame if available
        if (result == 0 && input.pos >= input.size) {
            // Try reading more to see if there's another frame
            size_t bytesRead = fread(inBuf, 1, ZSTD_DStreamInSize(), file);
            if (bytesRead > 0) {
                // Re-initialize for next frame
                ZSTD_initDStream(dstream);
                input.size = bytesRead;
                input.pos = 0;
                needMoreInput = true;
            } else {
                // No more data
                if (totalDecompressed > 0) {
                    break;  // We have data to return
                }
                reachedEOF = true;
            }
        }
    }

    if (totalDecompressed <= 0) // ERROR or EOF
        return EOF;

    // reset buffer pointers
    setg( buffer + (4 - n_putback),   // beginning of putback area
          buffer + 4,                 // read position
          buffer + 4 + totalDecompressed);  // end of buffer

    // return next character
    return * reinterpret_cast<unsigned char *>( gptr());
}

int zstdstreambuf::flush_buffer() {
    // Not used for input-only stream, but implemented for completeness
    int w = pptr() - pbase();
    pbump( -w);
    return w;
}

int zstdstreambuf::overflow( int c) { // used for output buffer only
    // Not implemented - FAME only needs input
    return EOF;
}

int zstdstreambuf::sync() {
    // For input stream, just return success
    if ( pptr() && pptr() > pbase()) {
        if ( flush_buffer() == EOF)
            return -1;
    }
    return 0;
}

// --------------------------------------
// class zstdstreambase:
// --------------------------------------

zstdstreambase::zstdstreambase( const char* name, int mode) {
    init( &buf);
    open( name, mode);
}

zstdstreambase::~zstdstreambase() {
    buf.close();
}

void zstdstreambase::open( const char* name, int open_mode) {
    if ( ! buf.open( name, open_mode))
        clear( rdstate() | std::ios::badbit);
}

void zstdstreambase::close() {
    if ( buf.is_open())
        if ( ! buf.close())
            clear( rdstate() | std::ios::badbit);
}

#ifdef ZSTDSTREAM_NAMESPACE
} // namespace ZSTDSTREAM_NAMESPACE
#endif

// ============================================================================
// EOF //