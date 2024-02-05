#ifndef TFACTORS_H
#define TFACTORS_H

struct tf_t {
    amrex::Real T9;
    amrex::Real T9i;
    amrex::Real T943i;
    amrex::Real T923i;
    amrex::Real T913i;
    amrex::Real T913;
    amrex::Real T923;
    amrex::Real T953;
    amrex::Real lnT9;
};

AMREX_GPU_HOST_DEVICE AMREX_INLINE
tf_t evaluate_tfactors(const amrex::Real T)
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
