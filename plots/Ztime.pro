PRO Ztime, opt 

window, xsize=700, ysize=700
TvLCT, 255, 0, 0, 2 
TvLCT, 128, 128, 128, 3 
Device, decomposed=0
!P.charsize=1.2

restore, 'halo_template.sav'

set_name = 'set71' 
; set_name = 'set125' 
FirstStarTime_data = read_ascii(strtrim(set_name)+'/halos_Firststartimecold.out', $
                                template=stars_template)
FirstStarTime = FirstStarTime_data.field001[1:*,490] ; yr 
HaloMass_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
HaloMass = HaloMass_data.field001[1:*,498];*1.0e10 ; yr 
ZcrTime_data = read_ascii(strtrim(set_name)+'/halos_Zcrtimecold.out', template=stars_template)
ZcrTime = ZcrTime_data.field001[1:*,490] ; yr 
EnrichmentTime = ZcrTime - FirstStarTime ; yr 

; set_name = 'set59' 
set_name = 'set165' 
FirstStarTime_data = read_ascii(strtrim(set_name)+'/halos_Firststartimecold.out', $
                                template=stars_template)
FirstStarTime2 = FirstStarTime_data.field001[1:*,90] ; yr 
HaloMass_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
HaloMass2 = HaloMass_data.field001[1:*,99]*1.0e10 ; yr 
ZcrTime_data = read_ascii(strtrim(set_name)+'/halos_Zcrtimecold.out', template=stars_template)
ZcrTime2 = ZcrTime_data.field001[1:*,90] ; yr 
EnrichmentTime2 = ZcrTime2 - FirstStarTime2 ; yr 

; set_name = 'set63' 
set_name = 'set171' 
FirstStarTime_data = read_ascii(strtrim(set_name)+'/halos_Firststartimecold.out', $
                                template=stars_template)
FirstStarTime3 = FirstStarTime_data.field001[1:*,90] ; yr 
HaloMass_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
HaloMass3 = HaloMass_data.field001[1:*,99]*1.0e10 ; yr 
ZcrTime_data = read_ascii(strtrim(set_name)+'/halos_Zcrtimecold.out', template=stars_template)
ZcrTime3 = ZcrTime_data.field001[1:*,90] ; yr 
EnrichmentTime3 = ZcrTime3 - FirstStarTime3 ; yr 

set_name = 'set70' 
; set_name = 'set126' 
FirstStarTime_data = read_ascii(strtrim(set_name)+'/halos_Firststartimecold.out', $
                                template=stars_template)
FirstStarTime4 = FirstStarTime_data.field001[1:*,490] ; yr 
HaloMass_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
HaloMass4 = HaloMass_data.field001[1:*,498]*1.0e10 ; yr 
ZcrTime_data = read_ascii(strtrim(set_name)+'/halos_Zcrtimecold.out', template=stars_template)
ZcrTime4 = ZcrTime_data.field001[1:*,490] ; yr 
EnrichmentTime4 = ZcrTime4 - FirstStarTime4 ; yr 

case opt of
   1: plot, FirstStarTime, /ylog
   2: plot, HaloMass, FirstStarTime, /ylog, /xlog, xtitle='!6halo mass', $
            ytitle='time of first star formation'
   3: begin 
      plot, HaloMass, FirstStarTime, /ylog, /xlog, xtitle='!6halo mass', $
            ytitle='time of first star formation'
      oplot, HaloMass, ZcrTime 
   end
   4: plot, EnrichmentTime, /ylog 
   5: plot, HaloMass, EnrichmentTime, /ylog, /xlog, yrange=[1.0e5,1.0e10], $
            xtitle='halo mass at z=0 (M!D!9n!X!N)', $
            ytitle='self-enrichment time scale (yr)'
   6: begin 
      plot, HaloMass, EnrichmentTime, /ylog, /xlog, yrange=[1.0e5,1.0e10], $
            xtitle='!6halo mass at z=0 (M!D!9n!X!N)', $
            ytitle='!6self-enrichment time scale (yr)', xrange=[1.0e9,1.0e12]
      foo = EnrichmentTime2[30:99]
      bar = rebin(foo, 280)
      ; oplot, HaloMass2, EnrichmentTime2, linestyle=2
      bar2 = bar[0:270]
      bar2 = smooth(bar2, 3)
      oplot, HaloMass2, bar2, linestyle=2
      oplot, HaloMass3, EnrichmentTime3, linestyle=3 
      legend, ['t!Ddelay!N=0 (fiducial)', 't!Ddelay!N=10!E8!Nyr', 't!Ddelay!N=10!E9!Nyr'], $
              linestyle=[0,2,3], charsize=1.1, /bottom
   end
   7: begin 
      set_plot, 'ps'
      device, filename='ztime.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0
      EnrichmentTime = EnrichmentTime4 + 3.0e6 
      EnrichmentTime = smooth(EnrichmentTime, 3) 
      plot, HaloMass, EnrichmentTime, /ylog, /xlog, yrange=[1.0e5,1.0e10], $
            xtitle='!6halo mass at z=0 (M!D!9n!X!N)', thick=3, $
            ytitle='!6self-enrichment time scale (yr)', xrange=[1.0e9,1.0e12]
      EnrichmentTime4 = smooth(EnrichmentTime4, 4) 
      oplot, HaloMass4, EnrichmentTime4, color=3, thick=3
      foo = EnrichmentTime2[30:99]
      bar = rebin(foo, 280)
      bar2 = bar[0:270]
      bar2 = smooth(bar2, 3)
      oplot, HaloMass2, bar2, linestyle=2, thick=3
      ; oplot, HaloMass2, EnrichmentTime2, linestyle=2, thick=3
      EnrichmentTime3 = smooth(EnrichmentTime3, 3) 
      oplot, HaloMass3, EnrichmentTime3, linestyle=3, thick=3 
      legend, ['t!Ddelay!N=0 (fiducial; 1-100 M!D!9n!X!N IMF)', $
               't!Ddelay!N=0 (fiducial; 100-260 M!D!9n!X!N IMF)', 't!Ddelay!N=10!E8!Nyr', $
               't!Ddelay!N=10!E9!Nyr'], $
              linestyle=[0,0,2,3], thick=[3,3,3,3], color=[-1,3,-1,-1], charsize=1.1, /bottom

      device, /close_file
      set_plot, 'X'
   end
   8: begin 
      readcol, '/home/girish/reion-eq/age.dat', h0t0, t0sec, t0yr 
      readcol, '/home/girish/reion-eq/z.dat', htt, zoft 
      local_hubble_0 = 1.023e-10*0.719
      zarr = zoft 
      tarr = t0yr(0)-htt/local_hubble_0
      FirstStarRedshift = interpol(zarr,tarr,FirstStarTime)
      len = size(zarr,/n_elements)

      readcol, 'deltat.chk', z, dz, dt 
      FirstStar_dt = interpol(dt,z,FirstStarRedshift)
      plot, HaloMass, EnrichmentTime, /ylog, /xlog, yrange=[1.0e5,1.0e10], $
            xtitle='!6halo mass at z=0 (M!D!9n!X!N)', $
            ytitle='!6self-enrichment time scale (yr)', xrange=[1.0e9,1.0e12]
      ;; plot, HaloMass, FirstStarTime, /ylog, /xlog, xtitle='!6halo mass', $
      ;;       ytitle='time of first star formation'
      oplot, HaloMass, FirstStar_dt, color=2
   end
   9: begin

      set_plot, 'ps'
      device, filename='ztime.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0

      plot, HaloMass, EnrichmentTime, /ylog, /xlog, yrange=[1.0e5,1.0e10], $
            xtitle='!6halo mass at z=0 (M!D!9n!X!N)', thick=3, $
            ytitle='!6self-enrichment time scale (yr)', xrange=[1.0e9,1.0e12]
      oplot, HaloMass2, EnrichmentTime2, linestyle=2, thick=3
      oplot, HaloMass3, EnrichmentTime3, linestyle=3, thick=3

      EnrichmentTime4 = smooth(EnrichmentTime4, 4) 
      oplot, HaloMass4, EnrichmentTime4, color=2, thick=3
      legend, ['t!Ddelay!N=0 (fiducial; low mass IMF)', $
               't!Ddelay!N=0 (fiducial; high mass IMF)', 't!Ddelay!N=10!E11!Nyr', $
               't!Ddelay!N=10!E12!Nyr'], $
              linestyle=[0,0,2,3], thick=[3,3,3,3], color=[-1,2,-1,-1], charsize=1.1, /bottom
      device, /close_file
      set_plot, 'X'

   end
endcase

END

