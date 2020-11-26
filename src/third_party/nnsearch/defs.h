#ifndef GEOGRAM_API_DEFS
#define GEOGRAM_API_DEFS

/**
 * \file geogram/api/defs.h
 * \brief Basic definitions for the Geogram C API
 */

/*
 * Deactivate warnings about documentation
 * We do that, because CLANG's doxygen parser does not know
 * some doxygen commands that we use (retval, copydoc) and
 * generates many warnings for them...
 */
#if defined(__clang__)
#pragma clang diagnostic ignored "-Wunknown-pragmas"
#pragma clang diagnostic ignored "-Wdocumentation-unknown-command"
#endif

/**
 * \brief Linkage declaration for geogram symbols.
 */

#if defined(GEO_DYNAMIC_LIBS)
   #if defined(_MSC_VER)
      #define GEO_IMPORT __declspec(dllimport)
      #define GEO_EXPORT __declspec(dllexport)
   #elif defined(__GNUC__)
      #define GEO_IMPORT
      #define GEO_EXPORT __attribute__ ((visibility("default")))
   #else
      #define GEO_IMPORT
      #define GEO_EXPORT
   #endif
#else
   #define GEO_IMPORT
   #define GEO_EXPORT
#endif

#ifdef geogram_EXPORTS
#define GEOGRAM_API GEO_EXPORT
#else
#define GEOGRAM_API GEO_IMPORT
#endif


/**
 * \brief A place-holder linkage declaration to indicate
 *  that the symbol should not be exported by Windows DLLs.
 * \details For instance, classes that inherit templates from
 *  the STL should not be exported, else it generates multiply
 *  defined symbols.
 */
#define NO_GEOGRAM_API

/**
 * \brief Opaque identifier of a mesh.
 * \details Used by the C API.
 */
typedef int GeoMesh;

/**
 * \brief Represents dimension (e.g. 3 for 3d, 4 for 4d ...).
 * \details Used by the C API.
 */
typedef unsigned char geo_coord_index_t;

/*
 * If GARGANTUA is defined, then geogram is compiled
 * with 64 bit indices.
 */
#ifdef GARGANTUA

#include <stdint.h>

/**
 * \brief Represents indices.
 * \details Used by the C API.
 */
typedef uint64_t geo_index_t;

/**
 * \brief Represents possibly negative indices.
 * \details Used by the C API.
 */
typedef int64_t geo_signed_index_t;

#else

/**
 * \brief Represents indices.
 * \details Used by the C API.
 */
typedef unsigned long geo_index_t;

/**
 * \brief Represents possibly negative indices.
 * \details Used by the C API.
 */
typedef int geo_signed_index_t;

#endif

/**
 * \brief Represents floating-point coordinates.
 * \details Used by the C API.
 */
typedef double geo_coord_t;

/**
 * \brief Represents truth values.
 * \details Used by the C API.
 */
typedef int geo_boolean;

/**
 * \brief Thruth values (geo_boolean).
 * \details Used by the C API.
 */
enum {
    GEO_FALSE = 0,
    GEO_TRUE = 1
};

#endif
