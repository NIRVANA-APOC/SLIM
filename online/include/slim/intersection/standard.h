#ifndef SLIM_INTERSECTION_STANDARD_H
#define SLIM_INTERSECTION_STANDARD_H

#include <cstdint>
#include <cstddef>

namespace slim
{

    inline size_t intersectStandard(const uint32_t *a, size_t a_len,
                                    const uint32_t *b, size_t b_len,
                                    uint32_t *result)
    {
        if (a_len == 0 || b_len == 0)
        {
            return 0;
        }

        size_t i = 0;
        size_t j = 0;
        size_t k = 0;

        while (i < a_len && j < b_len)
        {
            if (a[i] == b[j])
            {
                result[k++] = a[i];
                i++;
                j++;
            }
            else if (a[i] < b[j])
            {
                i++;
            }
            else
            {
                j++;
            }
        }

        return k;
    }

}

#endif
