#ifndef SLIM_INTERSECTION_AVX2_H
#define SLIM_INTERSECTION_AVX2_H

#include <cstdint>
#include <cstddef>
#include <algorithm>

#ifdef __AVX2__
#include <immintrin.h>
#endif

#include <slim/intersection/standard.h>

namespace slim
{

    inline size_t intersectAVX2(const uint32_t *a, size_t a_len,
                                const uint32_t *b, size_t b_len,
                                uint32_t *result)
    {
#ifdef __AVX2__
        if (a_len == 0 || b_len == 0)
        {
            return 0;
        }

        size_t lc = a_len;
        size_t rc = b_len;
        const uint32_t *larray = a;
        const uint32_t *rarray = b;

        if (lc > rc)
        {
            std::swap(larray, rarray);
            std::swap(lc, rc);
        }

        size_t li = 0;
        size_t ri = 0;
        size_t k = 0;

        auto size_ratio = static_cast<double>(rc) / lc;

        if (size_ratio > 2.0 && li < lc && ri + 7 < rc)
        {
            __m256i u_elements = _mm256_set1_epi32(larray[li]);
            __m256i v_elements = _mm256_loadu_si256(
                reinterpret_cast<const __m256i *>(rarray + ri));

            while (true)
            {
                __m256i mask = _mm256_cmpeq_epi32(u_elements, v_elements);
                int real_mask = _mm256_movemask_epi8(mask);

                if (real_mask != 0)
                {

                    for (int i = 0; i < 8 && ri + i < rc; ++i)
                    {
                        if (larray[li] == rarray[ri + i])
                        {
                            result[k++] = larray[li];
                            break;
                        }
                    }
                }

                if (larray[li] > rarray[ri + 7])
                {
                    ri += 8;
                    if (ri + 7 >= rc)
                        break;
                    v_elements = _mm256_loadu_si256(
                        reinterpret_cast<const __m256i *>(rarray + ri));
                }
                else
                {
                    li++;
                    if (li >= lc)
                        break;
                    u_elements = _mm256_set1_epi32(larray[li]);
                }
            }
        }

        while (li < lc && ri < rc)
        {
            if (larray[li] == rarray[ri])
            {
                result[k++] = larray[li];
                li++;
                ri++;
            }
            else if (larray[li] < rarray[ri])
            {
                li++;
            }
            else
            {
                ri++;
            }
        }

        return k;
#else

        return intersectStandard(a, a_len, b, b_len, result);
#endif
    }

}

#endif
