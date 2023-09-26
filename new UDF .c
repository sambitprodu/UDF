
# include "udf.h"
real source_mass_P1           /* p1 is liquid*/
real source_mass_P2         /* p2 is solid*/
real MINIMAL_DIA = 0.00001;							
real MINIMAL_VOF = 0.01;

real t_mix;
real t_P1;
real p2;
real t_liquid;
real t_solid;
real t_ms_P1;
real t_ms_P2;
real d;


/* defining the nucleation rate, this will be related to mass transfer further, LUDWIK EQN [21] */
DEFINE_SOURCE(nucleation, c, t_solid, dS, eqn)
{
	real T_CURRENT;
	real T_f;
	real T_prev;
	real liq_frac;
	real solid_frac;
	real n_max;
	double delta_H;
	double d_delta_T;
	double d_ms_P1;
	double d_n;
	double N1;
	int timestep;
	real CL_curr;
	
	
	t_mix = THREAD_SUPER_THREAD (t_solid);
	t_liquid = THREAD_SUB_THREAD (t_mix,0);
	T_CURRENT = C_T(c,t_liquid);
	T_prev = C_T_M1(c,t_liquid);
	liq_frac = C_VOF (C,t_liquid);
	CL_curr = C_UDSI(C,t_liquid,2);
	ds[eqn] =0.0;
	if(T_CURRENT>= 868 || liq_frac < 0.1)
	{
	N = 0.00;
	}
	else
	{
		T_f = 867;
		m = -290;
		n_max = 1*pow(10,14);
		delta_T_N = 20;
		delta_T_sigma = 5;
		delta_T = T_f+(m*CL_curr)-T_CURRENT;
		Time = (CURRENT_TIME - PREVIOUS_TIME);
		d_delta_T = (double)(T_prev-T_CURRENT)/Time;
		d_n = (double)((n_max/(sqrt(2.*3.14)*delta_T_sigma))*exp((-0.5*0.3)*pow(((delta_T-delta_T_N)/delta_T_sigma),2.)));
		timestep = RP_Get_Integer("time-step");
		N1 = (double)fabs(d_n*d_delta_T);
		}
		return N1;
		}

DEFINE_PROPERTY (Grain_Diameter_solid,c,t_solid)
{
	double Fs;
	double n;
	double Fl;
	real d;

	
	THREAD *t_mix;
	THREAD *t_liquid;
	t_mix = THREAD_SUPER_THREAD(t_solid);
	t_mix= THREAD_SUB_THREAD(t_mix,1);
	n = C_UDSI(c,t_solid);
	liq_frac = C_VOF (C,t_liquid);
	solid_frac = C_VOF (C,t_solid);
	if (liq_frac>= 0.5)
	{
		d =(double)pow(((24*soild_frac)/(n*4*3.14)),0.33);
	}
	else
	{
		d = 0.00001;
	}
	return d;
}
		
/* mass transfer rate FOR p1, LUDWIK EQN [18] */
DEFINE_ADJUST (mass_transfer_rate,c_ms_P1,t_ms_P1, d_ms_P1, eqn)
{
	real ro_p1;
	real ro_P2;
	real vof_P1;
	real vof_P2;
	real m;
	real G_alpha;
	real G;
	
	
	Thread *t_mix;
	Thread * t_liquid, *t_solid;
	real source_mass_P1;
	thread_loop_c (t_ms_P1, d_ms_P1);
	{
		t_P1=THREAD_SUB_THREAD(t_ms_P1, 0);
		t_P2=THREAD_SUB_THREAD(t_ms_P1, 1);
		if (FLUID_THREAD_P(t_liquid))
		{
			begin_c_loop(c_ms_P1, t_ms_P1);
			{
				ro_p1 = C_R(c_ms_P1, t_P1);
				ro_P2 = C_R(c_ms_P1, t_P2);
				vof_P1 = C_VOF (C_R(c_ms_P1, t_P1);
				vof_P2 = C_VOF (C_R(c_ms_P1, t_P2);
				T_f = 867;
				m = -290;
				G_alpha = 10;
				CL_curr = C_UDSI(C,t_liquid,2);
				G = G_alpha * (((T_CURRENT-T_f)/m)-CL_curr);
				if(C_VOF(c_ms_P1,t_liquid)<0.1)
				{
				source_mass_P1=0.0;
				}
				else
				{
					source_mass_P1 = G * 3.14 *d * d (c_ms_P1, t_solid)*(1-C_VOF(c_ms_P1,t_solid));
				}
				source_mass_P1 = C_UDMI ((c_ms_P1, t_ms_p1,0));
			}
			end_c_loop (c_ms_P1, t_ms_P1)
		}
	}
}
DEFINE_SOURCE (mass_P1, c_ms_P1, t_ms_P1, ds_ms_P1, eqn_ms_P1)
{
	THREAD *mix_t_ms_P2, *t_ms_P2;
	mix_t_ms_P1 = THREAD_SUPER_THREAD (t_ms_P1);
	t_ms_P2 = THREAD_SUB_THREAD(mix_t_ms_P1,1);
	ro_p1 = C_R(c_ms_P1, t_P1);
	ro_P2 = C_R(c_ms_P1, t_P2);
	vof_P1 = C_VOF (C_R(c_ms_P1, t_P1);
	vof_P2 = C_VOF (C_R(c_ms_P1, t_P2);
	timestep = RP_Get_Integer("time-step");
	ds_ms_P1[eqn_ms_P1]= 0.0;
	if(fabs(_UDMI(c_ms_P1,mix_t_ms_P1,3))0.0)
	return 0.0;
	if(vof_p1 <= MINIMAL_VOF)
	{
		source_mass_P1 = 0.0;
	}
	else
	{
		source_mass_P1=-1.0*( C_UDMI(c_ms_p1, mix_t_ms_p1,3);  /* if liquid phase is less than minimum volume then no phase transfer will be taking place*/
	
	}
	return source_mass_P1;
}

DEFINE_SOURCE (mass_P2, c_ms_P2, t_ms_P2, ds_ms_P2, eqn_ms_P2)
{
	THREAD *mix_t_ms_P2, *t_ms_P2;
	mix_t_ms_P1 = THREAD_SUPER_THREAD (t_ms_P1);
	t_ms_P2 = THREAD_SUB_THREAD(mix_t_ms_P1,1);
	ro_p1 = C_R(c_ms_P1, t_P1);
	ro_P2 = C_R(c_ms_P1, t_P2);
	vof_P1 = C_VOF (C_R(c_ms_P1, t_P1);
	vof_P2 = C_VOF (C_R(c_ms_P1, t_P2);
	timestep = RP_Get_Integer("time-step");
	ds_ms_P2[eqn_ms_P2]= 0.0;
	if(fabs(_UDMI(c_ms_P2,mix_t_ms_P2,4))0.0)
	return 0.0;
	if(vof_p2 <= MINIMAL_VOF)
	{
		source_mass_P1 = 0.0;
	}
	else
	{
		source_mass_P2=-1.0*( C_UDMI(c_ms_p2, mix_t_ms_p2,3);  /* if liquid phase is less than minimum volume then no phase transfer will be taking place*/
	
	}
	return source_mass_P2;
}

DEFINE_SOURCE(momentum_ls_P1_X,c,t_liquid,ds,eqn)
{
	u_l_x = C_U(c,t_liquid);
	u_s_x = C_U(c,t_liquid); 
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	ds[eqn] = 0;
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		u_l = u_l_x;
	}
	else
	{
		u_l = u_s_x;
	}
	U_ls_p_x = u_l*source_mass_P1;
	return U_ls_p_x:
}

DEFINE_SOURCE(momentum_ls_P1_Y,c,t_liquid,ds,eqn)
{
	u_l_y = C_U(c,t_liquid);
	u_s_y = C_U(c,t_liquid); 
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	ds[eqn] = 0;
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		u_l = u_l_y;
	}
	else
	{
		u_l = u_s_y;
	}
	U_ls_p_y =u_l*source_mass_P1;
	return U_ls_p_y:
}

DEFINE_SOURCE(momentum_ls_P1_Z,c,t_liquid,ds,eqn)
{
	u_l_z = C_U(c,t_liquid);
	u_s_z = C_U(c,t_liquid); 
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	ds[eqn] = 0;
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		u_l = u_l_z;
	}
	else
	{
		u_l = u_s_z;
	}
	U_ls_p_z =u_l*source_mass_P1;
	return U_ls_p_z:
}
	
DEFINE_SOURCE( Energy_source_ls,c,t_liquid,ds,eqn)
{
	t_mix = THREAD_SUPER_THREAD(t_liquid);
	t_solid = THREAD_SUB_THREAD(t_mix,1);
	ds[eqn]=0;
	liq_frac = C_VOF(c,t_liquid);
	solid_frac = C_VOF(c,t_solid);
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		h = H_l;
	}
	else
	{
	h = H_s
	}
	Q_ls = h* source_mass_P1;
	return Q_ls;
}

	DEFINE_SOURCE(momentum_sl_p2_X,c,t_solid, ds, eqn)
{
	u_l_x = C_U(c,t_liquid);
	u_s_x = C_U(c,t_liquid); 
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	ds[eqn] = 0;
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		u_l = u_l_x;
	}
	else
	{
		u_l = u_s_x;
	}
	
	U_sl_p_x = u_l*source_mass_P2;

	return U_sl_p_x;
}

DEFINE_SOURCE(momentum_sl_p2_Y,c,t_solid, ds, eqn)
{
	u_l_y = C_U(c,t_liquid);
	u_s_y = C_U(c,t_liquid); 
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	ds[eqn] = 0;
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		u_l = u_l_y;
	}
	else
	{
		u_l = u_s_y;
	}
	
	U_sl_p_y = u_l*source_mass_P2;

	return U_sl_p_y;
}

DEFINE_SOURCE(momentum_sl_p2_z,c,t_solid, ds, eqn)
{
	u_l_z = C_U(c,t_liquid);
	u_s_z = C_U(c,t_liquid); 
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	ds[eqn] = 0;
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		u_l = u_l_z;
	}
	else
	{
		u_l = u_s_z;
	}
	
	U_sl_p_z = u_l*source_mass_P2;

	return U_sl_p_z;
}
	
DEFINE_SOURCE( Energy_source_sl,c,t_solid,ds,eqn)        
{
	t_mix = THREAD_SUPER_THREAD(t_liquid);
	t_solid = THREAD_SUB_THREAD(t_mix,1);
	ds[eqn]=0;
	liq_frac = C_VOF(c,t_liquid);
	solid_frac = C_VOF(c,t_solid);
	H_l = C_H(c,t_liquid);
	H_s = C_H(c,t_solid);
	delta_H = H_l-H_s;
	if(delta_H>0)
	{
		h = H_l;
	}
	else
	{
	h = H_s
	}
	Q_sl= h* source_mass_P2;
	return Q_sl;
}

	
	
	
	
	
	
	
				
			
	
	

 
