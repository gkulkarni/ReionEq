PRO reion2, opt 

  case opt of

     1: begin 
     
     window, xsize=1000, ysize=1000
     Device, decomposed=0
     TvLCT, 255, 0, 0, 2 
     TvLCT, 0, 127, 255, 3
     TvLCT, 255, 255, 0, 4 
     !P.charsize = 2

     readcol, 'set161/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
              lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
  
     readcol, 'set201/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
              lmfp2, r, igmdcrit, nphdot2, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

     readcol, 'set202/reion.out', z3, q, tau, gammapi3, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
              lmfp2, r, igmdcrit, nphdot2, temphva, fv, igmdcrit, fm, $
              format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
  
     plot, z, gammapi, /xlog, /ylog, xrange=[0.1,100]
     oplot, z2, gammapi2, thick=4
     oplot, z3, gammapi3, thick=4, linestyle=5

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

     4: begin 
        
        window, xsize=1000, ysize=1000
        Device, decomposed=0
        TvLCT, 255, 0, 0, 2 
        TvLCT, 0, 127, 255, 3
        TvLCT, 255, 255, 0, 4 
        !P.charsize = 2

        readcol, 'set200/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
                 lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        
        readcol, 'set201/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph2, tempc, avtemp, x_ii, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva2, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        readcol, 'set202/reion.out', z3, q, tau, gammapi3, gpi_p2, gpi_p3, temph3, tempc, avtemp, x_ii, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva3, fv, igmdcrit, fm, $
                 format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        readcol, 'set171/reion.out', z4, q, tau, gammapi2, gpi_p2, gpi_p3, temph4, tempc, avtemp, x_ii, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva4, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        
        plot, z, temph, /xlog, /ylog, xrange=[0.1,100]
        oplot, z2, temph2, thick=4
        oplot, z3, temph3, thick=4, linestyle=5
        oplot, z4, temph4, color=2

     end

     5: begin 
        
        window, xsize=1000, ysize=1000
        Device, decomposed=0
        TvLCT, 255, 0, 0, 2 
        TvLCT, 0, 127, 255, 3
        TvLCT, 255, 255, 0, 4 
        !P.charsize = 2

        readcol, 'set200/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
                 lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        
        readcol, 'set201/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph2, tempc, avtemp, x_ii, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva2, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        readcol, 'set202/reion.out', z3, q, tau, gammapi3, gpi_p2, gpi_p3, temph3, tempc, avtemp, x_ii, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva3, fv, igmdcrit, fm, $
                 format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        readcol, 'set171/reion.out', z4, q, tau, gammapi2, gpi_p2, gpi_p3, temph4, tempc, avtemp, x_ii, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva4, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        
        plot, z, temphva, /xlog, /ylog, xrange=[0.1,100]
        oplot, z2, temphva2, thick=4
        oplot, z3, temphva3, thick=4, linestyle=5
        oplot, z4, temphva4, color=2

     end

     6: begin 
        
        window, xsize=1000, ysize=1000
        Device, decomposed=0
        TvLCT, 255, 0, 0, 2 
        TvLCT, 0, 127, 255, 3
        TvLCT, 255, 255, 0, 4 
        !P.charsize = 2

        readcol, 'set200/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
                 lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set201/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph2, tempc, avtemp, x_ii2, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva2, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set202/reion.out', z3, q, tau, gammapi3, gpi_p2, gpi_p3, temph3, tempc, avtemp, x_ii3, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva3, fv, igmdcrit, fm, $
                 format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set171/reion.out', z4, q, tau, gammapi2, gpi_p2, gpi_p3, temph4, tempc, avtemp, x_ii4, dnlldz, $
                 lmfp2, r, igmdcrit, nphdot2, temphva4, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        plot, z, x_ii, /xlog, /ylog, xrange=[0.1,100]
        oplot, z2, x_ii2, thick=4
        oplot, z3, x_ii3, thick=4, linestyle=5
        oplot, z4, x_ii4, color=2

     end

     7: begin 
        
        window, xsize=1000, ysize=1000
        Device, decomposed=0
        TvLCT, 255, 0, 0, 2 
        TvLCT, 0, 127, 255, 3
        TvLCT, 255, 255, 0, 4 
        !P.charsize = 2

        readcol, 'set200/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
                 lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set201/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph2, tempc, avtemp, x_ii2, dnlldz, $
                 lmfp2, r2, igmdcrit, nphdot2, temphva2, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set202/reion.out', z3, q, tau, gammapi3, gpi_p2, gpi_p3, temph3, tempc, avtemp, x_ii3, dnlldz, $
                 lmfp2, r3, igmdcrit, nphdot2, temphva3, fv, igmdcrit, fm, $
                 format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set170/reion.out', z4, q, tau, gammapi2, gpi_p2, gpi_p3, temph4, tempc, avtemp, x_ii4, dnlldz, $
                 lmfp2, r4, igmdcrit, nphdot2, temphva4, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        plot, z, r, /xlog, /ylog, xrange=[0.1,100]
        oplot, z2, r2, thick=4
        oplot, z3, r3, thick=4, linestyle=5
        oplot, z4, r4, color=2

     end

     8: begin 
        
        window, xsize=1000, ysize=1000
        Device, decomposed=0
        TvLCT, 255, 0, 0, 2 
        TvLCT, 0, 127, 255, 3
        TvLCT, 255, 255, 0, 4 
        !P.charsize = 2

        readcol, 'set200/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
                 lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set201/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph2, tempc2, avtemp, x_ii2, dnlldz, $
                 lmfp2, r2, igmdcrit, nphdot2, temphva2, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set202/reion.out', z3, q, tau, gammapi3, gpi_p2, gpi_p3, temph3, tempc3, avtemp, x_ii3, dnlldz, $
                 lmfp2, r3, igmdcrit, nphdot2, temphva3, fv, igmdcrit, fm, $
                 format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set172/reion.out', z4, q, tau, gammapi2, gpi_p2, gpi_p3, temph4, tempc4, avtemp, x_ii4, dnlldz, $
                 lmfp2, r4, igmdcrit, nphdot2, temphva4, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        plot, z, tempc, /xlog, /ylog, xrange=[0.1,100]
        oplot, z2, tempc2, thick=4
        oplot, z3, tempc3, thick=4, linestyle=5
        oplot, z4, tempc4, color=2

     end

     9: begin 
        
        window, xsize=1000, ysize=1000
        Device, decomposed=0
        TvLCT, 255, 0, 0, 2 
        TvLCT, 0, 127, 255, 3
        TvLCT, 255, 255, 0, 4 
        !P.charsize = 2

        readcol, 'set200/reion.out', z, q, tau, gammapi, gpi_p2, gpi_p3, temph, tempc, avtemp, x_ii, dnlldz, $
                 lmfp, r, igmdcrit, nphdot, temphva, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set201/reion.out', z2, q, tau, gammapi2, gpi_p2, gpi_p3, temph2, tempc2, avtemp, x_ii2, dnlldz, $
                 lmfp2, r2, igmdcrit2, nphdot2, temphva2, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set202/reion.out', z3, q, tau, gammapi3, gpi_p2, gpi_p3, temph3, tempc3, avtemp, x_ii3, dnlldz, $
                 lmfp2, r3, igmdcrit3, nphdot2, temphva3, fv, igmdcrit, fm, $
                 format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 
        readcol, 'set171/reion.out', z4, q, tau, gammapi2, gpi_p2, gpi_p3, temph4, tempc4, avtemp, x_ii4, dnlldz, $
                 lmfp2, r4, igmdcrit4, nphdot2, temphva4, fv, format='D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D,D', /silent 

        plot, z, igmdcrit, /xlog, /ylog, xrange=[0.1,100]
        oplot, z2, igmdcrit2, thick=4
        oplot, z3, igmdcrit3, thick=4, linestyle=5
        oplot, z4, igmdcrit4, color=2

     end


  endcase

END



