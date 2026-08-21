// pk_config.hpp
// Central project-wide configuration for PKD compression.
// Include this BEFORE any PKD headers.
#pragma once

// Define PK_IO_USE_ZLIB in the build system and link with zlib.
#ifdef PK_IO_USE_ZSTD
#  error "PK_IO_USE_ZSTD is unsupported; PKD files use zlib compression."
#endif
#ifndef PK_IO_USE_ZLIB
#  error "Define PK_IO_USE_ZLIB before including PKD headers."
#endif
