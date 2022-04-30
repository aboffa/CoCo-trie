#ifndef MAX_L_CONFIG_H_IN
#define MAX_L_CONFIG_H_IN

#ifdef BIG_ALPHABET  // ALPHABET_SIZE <= 96
    #define __MAX_L_THRS 18
    #define __MAX_L_THRS_EF 9
#else // ALPHABET_SIZE <= 27
    #define __MAX_L_THRS 25
    #define __MAX_L_THRS_EF 13
#endif

#endif// MAX_L_CONFIG_H_IN
