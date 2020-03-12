SUBDIRS =					\
	2_gbatch_04				\
	common					\
	4_xtrydis_04				\
	5_xprojtm_04				\
	6_xslope_09				\
	8_direct3d_04				\
	9_fireworks_04				\
	10_hiresspk_04				\
	11_gtune_04				\
	12a_random_shuf_05			\
	12_local_block_shift_04			\
	13_spkpat_04				\
	14_3djmpbg_04				\
	15_cyclesorter_04			\
	ogl_sb



curr_manver_pdf=Gravity_Manual_17-Oct-2017_rev_1.3.pdf
std_man_name=Gravity_Manual.pdf
curr_manver_doc=Gravity_Manual_17-Oct-2017_rev_1.3.doc

curr_flow_pdf=grav_cluster_flow_chart_v1.pdf
std_flow_name=GravityFlowchart.pdf

dist_doc_DATA = $(curr_manver_pdf) $(curr_flow_pdf) README	\
2-1-100steps-sig.jpg runsig_P.sh runsig_N.sh run_P.sh run_N.sh	\
param_P param_N 2004-7-29-re.gdt

EXTRA_DIST = $(curr_manver_doc) debian

install-exec-hook:
	-cd $(DESTDIR)$(bindir) && \
	$(LN_S) -f gbatch gsig

# symbolic link to manual
install-data-hook:
	-cd $(DESTDIR)$(docdir) && \
	$(LN_S) -f $(curr_manver_pdf) $(std_man_name)
	-cd $(DESTDIR)$(docdir) && \
	$(LN_S) -f $(curr_flow_pdf) $(std_flow_name)
