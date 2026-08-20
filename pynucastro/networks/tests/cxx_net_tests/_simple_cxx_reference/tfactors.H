#ifndef TFACTORS_H
#define TFACTORS_H

#include <cmath>
#include <amrex_bridge.H>

struct tf_t {
    Real T9;
    Real T9i;
    Real T943i;
    Real T923i;
    Real T913i;
    Real T913;
    Real T923;
    Real T953;
    Real lnT9;
};

inline
tf_t evaluate_tfactors(const Real T)
{

    tf_t tf;
    tf.T9 = T / 1.e9_rt;
    tf.T9i = 1.0_rt / tf.T9;
    tf.T913 = std::cbrt(tf.T9);
    tf.T913i = 1.0_rt / tf.T913;
    tf.T923i = tf.T913i * tf.T913i;
    tf.T943i = tf.T9i * tf.T913i;
    tf.T923 = tf.T913 * tf.T913;
    tf.T953 = tf.T9 * tf.T923;
    tf.lnT9 = std::log(tf.T9);

    return tf;
}

#endif