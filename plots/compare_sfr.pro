PRO compare_sfr, burst_input_file, nonburst_input_file, row, duty_cycle 

  ;; File: compare_sfr.pro 
  ;;  Cre: 2014-01-07 
  ;; 
  ;; This code was written to understand how the implementation of
  ;; starbursts in reion-eq/main.f90 works, and what values of the
  ;; starburst duty cycle are sensible.  It plots `burst_enhancement'
  ;; and `mstardot_halosc' arrays from reion-eq/main.f90, after
  ;; calling them `burst_data_at_selected_redshift' and
  ;; `nonburst_data_at_selected_redshift.' An example usage is
  ;; 
  ;; compare_sfr, "set174/halos_aux.out", "set174/halosc_sfrcontrib.out", 87, 0.1

  restore, 'halo_template.sav'

  burst_data = read_ascii(burst_input_file, template=stars_template)
  nonburst_data = read_ascii(nonburst_input_file, template=stars_template)

  lcolumn = 90 

  burst_data_at_selected_redshift = burst_data.field001(lcolumn:*, row) 
  nonburst_data_at_selected_redshift = nonburst_data.field001(lcolumn:*, row) 

  burst_data_at_selected_redshift = abs(burst_data_at_selected_redshift)
  burst_data_at_selected_redshift *= duty_cycle 

  plot, burst_data_at_selected_redshift, /ylog 
  oplot, nonburst_data_at_selected_redshift, linestyle=2

END


