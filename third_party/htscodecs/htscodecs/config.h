/* Minimal config.h for htscodecs - fqzcomp_qual only */
#ifndef CONFIG_H
#define CONFIG_H

#define HAVE_BUILTIN_PREFETCH 1

/* We don't use threading in htscodecs (rayon handles parallelism) */
/* #undef HAVE_PTHREAD */

#endif /* CONFIG_H */
