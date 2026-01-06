#ifndef SLIM_SET_INTERSECTION_H
#define SLIM_SET_INTERSECTION_H

#include <cstdint>
#include <cstddef>

#ifdef __AVX2__
#include <immintrin.h>
#include <x86intrin.h>
#define SLIM_HAS_AVX2 1
#else
#define SLIM_HAS_AVX2 0
#endif

namespace slim
{

#if SLIM_HAS_AVX2

    inline uint32_t BinarySearchForGallopingSearchAVX2(const uint32_t *array, uint32_t offset_beg, uint32_t offset_end, uint32_t val)
    {
        while (offset_end - offset_beg >= 16)
        {
            auto mid = static_cast<uint32_t>((static_cast<unsigned long>(offset_beg) + offset_end) / 2);
            _mm_prefetch((char *)&array[(static_cast<unsigned long>(mid + 1) + offset_end) / 2], _MM_HINT_T0);
            _mm_prefetch((char *)&array[(static_cast<unsigned long>(offset_beg) + mid) / 2], _MM_HINT_T0);
            if (array[mid] == val)
            {
                return mid;
            }
            else if (array[mid] < val)
            {
                offset_beg = mid + 1;
            }
            else
            {
                offset_end = mid;
            }
        }

        __m256i pivot_element = _mm256_set1_epi32(val);
        for (; offset_beg + 7 < offset_end; offset_beg += 8)
        {
            __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
            __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
            int mask = _mm256_movemask_epi8(cmp_res);
            if (mask != (int)0xffffffff)
            {
                return offset_beg + (_mm_popcnt_u32(mask) >> 2);
            }
        }
        if (offset_beg < offset_end)
        {
            auto left_size = offset_end - offset_beg;
            __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
            __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
            int mask = _mm256_movemask_epi8(cmp_res);
            int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
            mask &= cmp_mask;
            if (mask != cmp_mask)
            {
                return offset_beg + (_mm_popcnt_u32(mask) >> 2);
            }
        }
        return offset_end;
    }

    inline uint32_t GallopingSearchAVX2(const uint32_t *array, uint32_t offset_beg, uint32_t offset_end, uint32_t val)
    {
        if (array[offset_end - 1] < val)
        {
            return offset_end;
        }

        __m256i pivot_element = _mm256_set1_epi32(val);
        if (offset_end - offset_beg >= 8)
        {
            __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
            __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
            int mask = _mm256_movemask_epi8(cmp_res);
            if (mask != (int)0xffffffff)
            {
                return offset_beg + (_mm_popcnt_u32(mask) >> 2);
            }
        }
        else
        {
            auto left_size = offset_end - offset_beg;
            __m256i elements = _mm256_loadu_si256(reinterpret_cast<const __m256i *>(array + offset_beg));
            __m256i cmp_res = _mm256_cmpgt_epi32(pivot_element, elements);
            int mask = _mm256_movemask_epi8(cmp_res);
            int cmp_mask = 0xffffffff >> ((8 - left_size) << 2);
            mask &= cmp_mask;
            if (mask != cmp_mask)
            {
                return offset_beg + (_mm_popcnt_u32(mask) >> 2);
            }
        }

        auto jump_idx = 8u;
        while (true)
        {
            auto peek_idx = offset_beg + jump_idx;
            if (peek_idx >= offset_end)
            {
                return BinarySearchForGallopingSearchAVX2(array, (jump_idx >> 1) + offset_beg + 1, offset_end, val);
            }
            if (array[peek_idx] < val)
            {
                jump_idx <<= 1;
            }
            else
            {
                return array[peek_idx] == val ? peek_idx : BinarySearchForGallopingSearchAVX2(array, (jump_idx >> 1) + offset_beg + 1, peek_idx + 1, val);
            }
        }
    }

    inline size_t ComputeCNGallopingAVX2(const uint32_t *larray, uint32_t l_count,
                                         const uint32_t *rarray, uint32_t r_count,
                                         uint32_t *cn)
    {
        size_t cn_count = 0;

        if (l_count == 0 || r_count == 0)
            return 0;

        uint32_t lc = l_count;
        uint32_t rc = r_count;

        if (lc > rc)
        {
            auto tmp = larray;
            larray = rarray;
            rarray = tmp;

            uint32_t tmp_count = lc;
            lc = rc;
            rc = tmp_count;
        }

        uint32_t li = 0;
        uint32_t ri = 0;

        while (true)
        {
            while (larray[li] < rarray[ri])
            {
                li += 1;
                if (li >= lc)
                {
                    return cn_count;
                }
            }

            ri = GallopingSearchAVX2(rarray, ri, rc, larray[li]);
            if (ri >= rc)
            {
                return cn_count;
            }

            if (larray[li] == rarray[ri])
            {
                cn[cn_count++] = larray[li];
                li += 1;
                ri += 1;
                if (li >= lc || ri >= rc)
                {
                    return cn_count;
                }
            }
        }
    }

    inline size_t ComputeCNMergeBasedAVX2(const uint32_t *larray, uint32_t l_count,
                                          const uint32_t *rarray, uint32_t r_count,
                                          uint32_t *cn)
    {
        size_t cn_count = 0;

        if (l_count == 0 || r_count == 0)
            return 0;

        uint32_t lc = l_count;
        uint32_t rc = r_count;

        if (lc > rc)
        {
            auto tmp = larray;
            larray = rarray;
            rarray = tmp;

            uint32_t tmp_count = lc;
            lc = rc;
            rc = tmp_count;
        }

        uint32_t li = 0;
        uint32_t ri = 0;

        __m256i per_u_order = _mm256_set_epi32(1, 1, 1, 1, 0, 0, 0, 0);
        __m256i per_v_order = _mm256_set_epi32(3, 2, 1, 0, 3, 2, 1, 0);
        uint32_t *cur_back_ptr = cn;

        auto size_ratio = (rc) / (lc);
        if (size_ratio > 2)
        {
            if (li < lc && ri + 7 < rc)
            {
                __m256i u_elements = _mm256_set1_epi32(larray[li]);
                __m256i v_elements = _mm256_loadu_si256((__m256i *)(rarray + ri));

                while (true)
                {
                    __m256i mask = _mm256_cmpeq_epi32(u_elements, v_elements);
                    auto real_mask = _mm256_movemask_epi8(mask);
                    if (real_mask != 0)
                    {
                        *cur_back_ptr = larray[li];
                        cur_back_ptr += 1;
                    }
                    if (larray[li] > rarray[ri + 7])
                    {
                        ri += 8;
                        if (ri + 7 >= rc)
                        {
                            break;
                        }
                        v_elements = _mm256_loadu_si256((__m256i *)(rarray + ri));
                    }
                    else
                    {
                        li++;
                        if (li >= lc)
                        {
                            break;
                        }
                        u_elements = _mm256_set1_epi32(larray[li]);
                    }
                }
            }
        }
        else
        {
            if (li + 1 < lc && ri + 3 < rc)
            {
                __m256i u_elements = _mm256_loadu_si256((__m256i *)(larray + li));
                __m256i u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                __m256i v_elements = _mm256_loadu_si256((__m256i *)(rarray + ri));
                __m256i v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);

                while (true)
                {
                    __m256i mask = _mm256_cmpeq_epi32(u_elements_per, v_elements_per);
                    auto real_mask = _mm256_movemask_epi8(mask);
                    if ((real_mask << 16) != 0)
                    {
                        *cur_back_ptr = larray[li];
                        cur_back_ptr += 1;
                    }
                    if ((real_mask >> 16) != 0)
                    {
                        *cur_back_ptr = larray[li + 1];
                        cur_back_ptr += 1;
                    }

                    if (larray[li + 1] == rarray[ri + 3])
                    {
                        li += 2;
                        ri += 4;
                        if (li + 1 >= lc || ri + 3 >= rc)
                        {
                            break;
                        }
                        u_elements = _mm256_loadu_si256((__m256i *)(larray + li));
                        u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                        v_elements = _mm256_loadu_si256((__m256i *)(rarray + ri));
                        v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                    }
                    else if (larray[li + 1] > rarray[ri + 3])
                    {
                        ri += 4;
                        if (ri + 3 >= rc)
                        {
                            break;
                        }
                        v_elements = _mm256_loadu_si256((__m256i *)(rarray + ri));
                        v_elements_per = _mm256_permutevar8x32_epi32(v_elements, per_v_order);
                    }
                    else
                    {
                        li += 2;
                        if (li + 1 >= lc)
                        {
                            break;
                        }
                        u_elements = _mm256_loadu_si256((__m256i *)(larray + li));
                        u_elements_per = _mm256_permutevar8x32_epi32(u_elements, per_u_order);
                    }
                }
            }
        }

        cn_count = (size_t)(cur_back_ptr - cn);
        if (li < lc && ri < rc)
        {
            while (true)
            {
                while (larray[li] < rarray[ri])
                {
                    ++li;
                    if (li >= lc)
                    {
                        return cn_count;
                    }
                }
                while (larray[li] > rarray[ri])
                {
                    ++ri;
                    if (ri >= rc)
                    {
                        return cn_count;
                    }
                }
                if (larray[li] == rarray[ri])
                {
                    cn[cn_count++] = larray[li];
                    ++li;
                    ++ri;
                    if (li >= lc || ri >= rc)
                    {
                        return cn_count;
                    }
                }
            }
        }
        return cn_count;
    }

    inline size_t intersect_adaptive(
        const uint32_t *a, size_t a_size,
        const uint32_t *b, size_t b_size,
        uint32_t *out)
    {
        if (a_size == 0 || b_size == 0)
            return 0;

        if (a_size / 50 > b_size || b_size / 50 > a_size)
        {
            return ComputeCNGallopingAVX2(a, (uint32_t)a_size, b, (uint32_t)b_size, out);
        }
        else
        {
            return ComputeCNMergeBasedAVX2(a, (uint32_t)a_size, b, (uint32_t)b_size, out);
        }
    }

#else

    inline size_t intersect_scalar(
        const uint32_t *a, size_t a_size,
        const uint32_t *b, size_t b_size,
        uint32_t *out)
    {
        size_t i = 0, j = 0, k = 0;
        while (i < a_size && j < b_size)
        {
            if (a[i] < b[j])
            {
                ++i;
            }
            else if (a[i] > b[j])
            {
                ++j;
            }
            else
            {
                out[k++] = a[i];
                ++i;
                ++j;
            }
        }
        return k;
    }

    inline size_t intersect_adaptive(
        const uint32_t *a, size_t a_size,
        const uint32_t *b, size_t b_size,
        uint32_t *out)
    {
        return intersect_scalar(a, a_size, b, b_size, out);
    }

#endif

}

#endif
