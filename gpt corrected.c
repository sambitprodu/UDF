#include "udf.h"

real MINIMAL_DIA = 0.00001;
real MINIMAL_VOF = 0.01;
real t_mix;
real t_P1;
real t_P2;
real t_liquid;
real t_solid;
real t_ms_P1;
real t_ms_P2;
real d;
real source_mass_P1;
real source_mass_P2;

/* Define properties */
DEFINE_PROPERTY(Grain_Diameter_solid, c, t_solid)
{
    real Fs;
    real n;
    real Fl;

    real d;
    THREAD *t_mix;
    THREAD *t_liquid;

    t_mix = THREAD_SUPER_THREAD(t_solid);
    t_liquid = THREAD_SUB_THREAD(t_mix, 0);

    n = C_UDSI(c, t_solid);
    real liq_frac = C_VOF(c, t_liquid);
    real solid_frac = C_VOF(c, t_solid);

    if (liq_frac >= 0.5)
    {
        d = pow((24 * solid_frac) / (n * 4 * 3.14), 0.333);
    }
    else
    {
        d = MINIMAL_DIA;
    }
    return d;
}

/* Mass transfer rate for P1 */
DEFINE_ADJUST(mass_transfer_rate, c_ms_P1, t_ms_P1, d_ms_P1, eqn)
{
    real ro_P1;
    real ro_P2;
    real vof_P1;
    real vof_P2;
    real G_alpha;
    real G;

    Thread *t_mix;
    Thread *t_liquid, *t_solid;

    thread_loop_c(t_ms_P1, d_ms_P1)
    {
        t_P1 = THREAD_SUB_THREAD(t_ms_P1, 0);
        t_P2 = THREAD_SUB_THREAD(t_ms_P1, 1);

        if (FLUID_THREAD_P(t_liquid))
        {
            begin_c_loop(c_ms_P1, t_ms_P1)
            {
                ro_P1 = C_R(c_ms_P1, t_P1);
                ro_P2 = C_R(c_ms_P1, t_P2);
                vof_P1 = C_VOF(c_ms_P1, t_P1);
                vof_P2 = C_VOF(c_ms_P1, t_P2);

                real T_f = 867.0;
                real m = -290.0;
                real CL_curr = C_UDSI(c_ms_P1, t_liquid, 2);

                G_alpha = 10.0;
                G = G_alpha * (((C_T(c_ms_P1, t_liquid) - T_f) / m) - CL_curr);

                if (C_VOF(c_ms_P1, t_liquid) < 0.1)
                {
                    source_mass_P1 = 0.0;
                }
                else
                {
                    source_mass_P1 = G * 3.14 * d * d * (1.0 - C_VOF(c_ms_P1, t_solid));
                }

                C_UDMI(c_ms_P1, t_ms_P1, 0) = source_mass_P1;
            }
            end_c_loop(c_ms_P1, t_ms_P1)
        }
    }
}

/* Define source terms for P1 and P2 */
DEFINE_SOURCE(mass_P1, c_ms_P1, t_ms_P1, ds_ms_P1, eqn_ms_P1)
{
    THREAD *mix_t_ms_P1, *t_ms_P2;
    real ro_P1, ro_P2;
    real vof_P1, vof_P2;
    real source_mass_P1;
    real timestep = RP_Get_Real("time-step");

    mix_t_ms_P1 = THREAD_SUPER_THREAD(t_ms_P1);
    t_ms_P2 = THREAD_SUB_THREAD(mix_t_ms_P1, 1);

    ro_P1 = C_R(c_ms_P1, t_P1);
    ro_P2 = C_R(c_ms_P1, t_P2);
    vof_P1 = C_VOF(c_ms_P1, t_P1);
    vof_P2 = C_VOF(c_ms_P1, t_P2);

    ds_ms_P1[eqn_ms_P1] = 0.0;

    if (fabs(C_UDMI(c_ms_P1, mix_t_ms_P1, 3)) > 0.0)
    {
        return 0.0;
    }

    if (vof_P1 <= MINIMAL_VOF)
    {
        source_mass_P1 = 0.0;
    }
    else
    {
        source_mass_P1 = -1.0 * C_UDMI(c_ms_P1, mix_t_ms_P1, 3); /* if liquid phase is less than minimum volume then no phase transfer will be taking place */
    }

    return source_mass_P1;
}

DEFINE_SOURCE(mass_P2, c_ms_P2, t_ms_P2, ds_ms_P2, eqn_ms_P2)
{
    THREAD *mix_t_ms_P2, *t_ms_P2;
    real ro_P1, ro_P2;
    real vof_P1, vof_P2;
    real source_mass_P2;
    real timestep = RP_Get_Real("time-step");

    mix_t_ms_P2 = THREAD_SUPER_THREAD(t_ms_P2);
    t_ms_P2 = THREAD_SUB_THREAD(mix_t_ms_P2, 1);

    ro_P1 = C_R(c_ms_P2, t_P1);
    ro_P2 = C_R(c_ms_P2, t_P2);
    vof_P1 = C_VOF(c_ms_P2, t_P1);
    vof_P2 = C_VOF(c_ms_P2, t_P2);

    ds_ms_P2[eqn_ms_P2] = 0.0;

    if (fabs(C_UDMI(c_ms_P2, mix_t_ms_P2, 4)) > 0.0)
    {
        return 0.0;
    }

    if (vof_P2 <= MINIMAL_VOF)
    {
        source_mass_P2 = 0.0;
    }
    else
    {
        source_mass_P2 = -1.0 * C_UDMI(c_ms_P2, mix_t_ms_P2, 4); /* if liquid phase is less than minimum volume then no phase transfer will be taking place */
    }

    return source_mass_P2;
}

/* Define momentum source terms */
DEFINE_SOURCE(momentum_ls_P1_X, c, t_liquid, ds, eqn)
{
    real u_l_x, u_s_x, H_l, H_s, delta_H, u_l, U_ls_p_x;

    u_l_x = C_U(c, t_liquid);
    u_s_x = C_U(c, t_liquid);
    H_l = C_H(c, t_liquid);
    H_s = C_H(c, t_solid);
    ds[eqn] = 0;
    delta_H = H_l - H_s;

    if (delta_H > 0)
    {
        u_l = u_l_x;
    }
    else
    {
        u_l = u_s_x;
    }

    U_ls_p_x = u_l * source_mass_P1;
    return U_ls_p_x;
}

/* Similar functions for Y and Z directions */

/* Energy source term for P1 */
DEFINE_SOURCE(Energy_source_ls, c, t_liquid, ds, eqn)
{
    real H_l, H_s, delta_H, h, Q_ls;

    t_mix = THREAD_SUPER_THREAD(t_liquid);
    t_solid = THREAD_SUB_THREAD(t_mix, 1);
    ds[eqn] = 0;
    real liq_frac = C_VOF(c, t_liquid);
    real solid_frac = C_VOF(c, t_solid);
    H_l = C_H(c, t_liquid);
    H_s = C_H(c, t_solid);
    delta_H = H_l - H_s;

    if (delta_H > 0)
    {
        h = H_l;
    }
    else
    {
        h = H_s;
    }

    Q_ls = h * source_mass_P1;
    return Q_ls;
}

/* Similar functions for momentum source terms for P2 */

/* Energy source term for P2 */

DEFINE_SOURCE(Energy_source_sl, c, t_solid, ds, eqn)
{
    real H_l, H_s, delta_H, h, Q_sl;

    t_mix = THREAD_SUPER_THREAD(t_solid);
    t_solid = THREAD_SUB_THREAD(t_mix, 1);
    ds[eqn] = 0;
    real liq_frac = C_VOF(c, t_liquid);
    real solid_frac = C_VOF(c, t_solid);
    H_l = C_H(c, t_liquid);
    H_s = C_H(c, t_solid);
    delta_H = H_l - H_s;

    if (delta_H > 0)
    {
        h = H_l;
    }
    else
    {
        h = H_s;
    }

    Q_sl = h * source_mass_P2;
    return Q_sl;
}
