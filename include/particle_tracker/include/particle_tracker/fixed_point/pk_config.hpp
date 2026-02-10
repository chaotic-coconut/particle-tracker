// pk_config.hpp
// Central project-wide configuration for PKD compression backend.
// Include this BEFORE any PKD headers.
#pragma once

// Choose exactly ONE codec in your build system instead of here:
//   -DPK_IO_USE_ZLIB   (link with -lz)
//   -DPK_IO_USE_ZSTD   (link with -lzstd)
//
// #define PK_IO_USE_ZLIB
// #define PK_IO_USE_ZSTD

#if !defined(PK_IO_USE_ZSTD) && !defined(PK_IO_USE_ZLIB)
#  error "Define either PK_IO_USE_ZSTD or PK_IO_USE_ZLIB before including PKD headers."
#endif
#if defined(PK_IO_USE_ZSTD) && defined(PK_IO_USE_ZLIB)
#  error "Define only ONE of PK_IO_USE_ZSTD or PK_IO_USE_ZLIB, not both."
#endif

