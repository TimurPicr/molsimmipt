# 3d Lennard-Jones melt

variable        rho index 0.85
variable        temp index 1.5
variable        sim_steps index 5000
log             log.diffusion.msd.T${temp}.rho${rho}
units		lj
atom_style	atomic

pair_style	lj/smooth 2.5 3.0

read_data       liquid.T${temp}.rho${rho}.data

neighbor	0.3 bin
neigh_modify	every 20 delay 0 check no

fix		1 all nve

timestep	0.002

reset_timestep  0

compute		msd all msd
variable        msd_xyz equal c_msd[4]
variable        step0 equal step
variable        step0 equal ${step0}
variable        dtime equal (step-${step0})*dt
variable        diffus equal v_msd_xyz/(6*(v_dtime+1e-12))
thermo_style	custom step time temp ke pe etotal press c_msd[4]
thermo          50

fix             print_msd all print 50 "${dtime} ${msd_xyz} ${diffus}" file msd.T${temp}.rho${rho}.log screen no
run		${sim_steps}
print           "Diffusivity estimate (MSD): ${diffus}" append msd.T${temp}.rho${rho}.log

write_data      liquid.T${temp}.rho${rho}.data
