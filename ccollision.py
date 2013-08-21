#!/usr/bin/env python

"""
ccollision.py

SciPy weave compiled collision testers.

"""


from scipy import weave
from scipy.weave import converters


def wdist(c_coa, f_coa, c_cob, f_cob, box):
    """Rewritten min_vect for weave speedup."""
    code = """
    int i, j, k;
    double p, q, r;
    double distance;
    double diff[3];
    double new_f_cob[3];
    int no_min_image;

    no_min_image = 0;

    for(i=0;i<3;i++)
    {
        p = f_coa[i];
        q = f_cob[i];
        diff[i] = p - q;
        if (diff[i] < -0.5)
        {
            new_f_cob[i] = q - 1.0;
        }
        else
        {
            if (diff[i] > 0.5)
            {
                new_f_cob[i] = q + 1.0;
            }
            else
            {
                new_f_cob[i] = q;
                no_min_image++;
            };
        };
    };

//    if (diff[0] != 0.0 || diff[1] != 0.0 || diff[2] != 0.0)
    if (no_min_image < 3)
    {
//        std::cout << diff[0] << diff[1] << diff[2] << std::endl;

        distance = 0;

        for(i=0;i<3;i++)
        {
            p = new_f_cob[0]*box(0, i) + new_f_cob[1]*box(1, i) + new_f_cob[2]*box(2, i);
            q = c_coa(i);
            distance = distance + pow(p - q, 2);
        };

    }
    else
    {
//        std::cout << diff[0] << std::endl;

        distance = 0;

        for(i=0;i<3;i++)
        {
            p = c_cob(i);
            q = c_coa(i);
            distance = distance + pow(p - q, 2);
        };
    };

    return_val = sqrt(distance);
    """
    return weave.inline(code, ['c_coa', 'f_coa', 'c_cob', 'f_cob', 'box'],
                        type_converters=converters.blitz,
                        support_code='#include <math.h>')
