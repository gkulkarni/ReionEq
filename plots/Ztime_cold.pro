PRO Ztime_cold, opt 

window, xsize=700, ysize=700
Device, decomposed=0
!P.charsize=1.2

restore, 'halo_template.sav'

set_name = 'set60' 

FirstStarTime_data = read_ascii(strtrim(set_name)+'/halos_Firststartimecold.out', $
                                template=stars_template)
FirstStarTime = FirstStarTime_data.field001[1:*,90] ; yr 
 
HaloMass_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
HaloMass = HaloMass_data.field001[1:*,99]*1.0e10 ; yr 

ZcrTime_data = read_ascii(strtrim(set_name)+'/halos_Zcrtimecold.out', template=stars_template)
ZcrTime = ZcrTime_data.field001[1:*,90] ; yr 

EnrichmentTime = ZcrTime - FirstStarTime ; yr 

set_name = 'set59' 

FirstStarTime_data = read_ascii(strtrim(set_name)+'/halos_Firststartimecold.out', $
                                template=stars_template)
FirstStarTime2 = FirstStarTime_data.field001[1:*,90] ; yr 
 
HaloMass_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
HaloMass2 = HaloMass_data.field001[1:*,99]*1.0e10 ; yr 

ZcrTime_data = read_ascii(strtrim(set_name)+'/halos_Zcrtimecold.out', template=stars_template)
ZcrTime2 = ZcrTime_data.field001[1:*,90] ; yr 

EnrichmentTime2 = ZcrTime2 - FirstStarTime2 ; yr 

set_name = 'set63' 

FirstStarTime_data = read_ascii(strtrim(set_name)+'/halos_Firststartimecold.out', $
                                template=stars_template)
FirstStarTime3 = FirstStarTime_data.field001[1:*,90] ; yr 
 
HaloMass_data = read_ascii(strtrim(set_name)+'/halos.out', template=stars_template)
HaloMass3 = HaloMass_data.field001[1:*,99]*1.0e10 ; yr 

ZcrTime_data = read_ascii(strtrim(set_name)+'/halos_Zcrtimecold.out', template=stars_template)
ZcrTime3 = ZcrTime_data.field001[1:*,90] ; yr 

EnrichmentTime3 = ZcrTime3 - FirstStarTime3 ; yr 

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
      oplot, HaloMass2, EnrichmentTime2, linestyle=2
      oplot, HaloMass3, EnrichmentTime3, linestyle=1 
   end
   7: begin 
      set_plot, 'ps'
      device, filename='ztime.ps', xsize=7.0, ysize=7.0, /inches, color=1, yoffset=1.0
      plot, HaloMass, EnrichmentTime, /ylog, /xlog, yrange=[1.0e5,1.0e10], $
            xtitle='!6halo mass at z=0 (M!D!9n!X!N)', $
            ytitle='!6self-enrichment time scale (yr)', xrange=[1.0e8,1.0e16]
      oplot, HaloMass2, EnrichmentTime2, linestyle=2
      oplot, HaloMass3, EnrichmentTime3, linestyle=1 
      device, /close_file
      set_plot, 'X'
end
endcase

END

