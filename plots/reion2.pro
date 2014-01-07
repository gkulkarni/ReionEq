PRO reion2, opt 

  case opt of

     1: begin 
     
     window, xsize=1500, ysize=500
     Device, decomposed=0
     TvLCT, 255, 0, 0, 2 
     TvLCT, 0, 127, 255, 3
     TvLCT, 255, 255, 0, 4 
     !P.charsize = 2

     readcol, 'set135/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
              lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
  
     readcol, 'set134/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
              lmfp2, r, igmdcrit, nphdot2, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
  
     plot, z, nphdot, /xlog, /ylog, xrange=[10,100]
     oplot, z2, nphdot2, color=2

     plot, z, lmfp, /xlog, /ylog, xrange=[10,100]
     oplot, z2, lmfp2, color=2

     plot, z, gammapi, /xlog, /ylog, xrange=[10,100]
     oplot, z2, gammapi2, color=2

  end

     2: begin 

     window, xsize=1000, ysize=1000
     Device, decomposed=0
     TvLCT, 255, 0, 0, 2 
     TvLCT, 0, 127, 255, 3
     TvLCT, 255, 255, 0, 4 
     !P.charsize = 2

     readcol, 'set135/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
              lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
  
     readcol, 'set134/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
              lmfp2, r, igmdcrit, nphdot2, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
     
     sample = gammapi[83:96] 
     sample = smooth(sample,7) 
     for i = 83, 96 do begin 
        gammapi[i] = sample[i-83] 
     endfor
     plot, z, gammapi, /xlog, /ylog, xrange=[1,100]
          
     ;; for i = 0, size(gammapi,/n_elements)-1 do begin 
     ;;    print, i, z[i], gammapi[i] 
     ;; endfor
     sample = gammapi2[83:99] 
     sample = smooth(sample,6) 
     for i = 83, 99 do begin 
        gammapi2[i] = sample[i-83] 
     endfor

     oplot, z2, gammapi2, color=2

  end

     3: begin 

        window, xsize=1000, ysize=1000
        Device, decomposed=0
        TvLCT, 255, 0, 0, 2 
        TvLCT, 0, 127, 255, 3
        TvLCT, 255, 255, 0, 4 
        !P.charsize = 2

        readcol, 'set135/sfr.out', z, sfr, sfr2, sfr3, lim, /silent
        readcol, 'set134/sfr.out', z2, sfr2, sfr22, sfr32, lim2, /silent 
        
        sample = sfr[86:96] 
        sample = smooth(sample,3) 
        for i = 86, 96 do begin 
           sfr[i] = sample[i-86] 
        endfor
        plot, z, sfr, /xlog, /ylog, xrange=[1,100], yrange=[1.0e-6,1]

        sample = sfr2[86:99] 
        sample = smooth(sample,3) 
        for i = 86, 99 do begin 
           sfr2[i] = sample[i-86] 
        endfor
        oplot, z2, sfr2, color=2
        
     end
     
  endcase

END



