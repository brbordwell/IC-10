pro run_em, ell = ell, mc = mc, combo = combo
;this procedure will run all the possible combinations of options for
;the isochrone fit
;/combo = does core+IR instead of sky

  common ic10, info
  difs = [.08, .1, .1] ;providing the cutoffs for contaminants (empirical)
  weight = .25 ;providing the weight of FG contamination (empirical)
  mem_weight = 4  ;providing a weight for the absolutely proven members (empirical)
  info = {center: [20/15.+17.3/240., 59+18/60.+13.6/3600.], rad:3/60., axrat: .61, posangl:-48., dif: difs, cont_weight: weight, mem_weight: mem_weight}
  ;[center/rad/posangl]= degrees
  ;[dif] = magnitudes
  ;[cont_weight, mem_weight] = unitless

  if n_elements(file_search('stellar*_pops.sav')) NE 4 then pop_main
  if n_elements(file_search('foreground*.sav')) NE 4 then fg_main
  ;if alternate structures of rings/skyrings/angles exist, this is
  ;...where it will be manually accounted for using the altrings option

  if keyword_set(combo) then begin
     if keyword_set(ell) then begin
        if keyword_set(mc) then begin
           iso_main, /ellipse, /masscut, /combo
           iso_main, /ours, /ellipse, /masscut, /combo
        endif else begin
           iso_main, /ellipse, /combo
           iso_main, /ours, /ellipse, /combo
        endelse
     endif else begin
        if keyword_set(mc) then begin
           iso_main, /masscut, /combo
           iso_main, /ours, /masscut, /combo
        endif else begin
           iso_main, /combo
           iso_main, /ours, /combo
        endelse
     endelse
  endif else print, "No sky this night"
end




pro iso_main, ellipse = ellipse, masscut = masscut, ours = ours, altrings = altrings, combo = combo
;+
;PURPOSE:

;INPUT:
;/ellipse = use the elliptical population cutoffs
;/masscut = only use isochrones points with a mass equivalent to that of
;a supergiant
;/ours = use the population cutoffs based on the location of our
;spectroscopic members
;/altrings = the user may provide a string with an alternate set 
;of rings/skyrings/angular cutoffs for the population divisions than the default

;OUTPUT:
;-
  common ic10

  if keyword_set(altrings) then bit = altrings else bit = ''
  if keyword_set(ellipse) then begin
     if keyword_set(ours) then bit += '_e_m' else bit += '_e'
  endif else begin
     if keyword_set(ours) then bit += '_m' else bit += ''
  endelse

  fil = file_search('stellar'+bit+'_pops.sav')  &  restore, fil
  fg_fil = file_search('foreground'+bit+'.sav')  &  restore, fg_fil
  ;restoring a structure (pops) with sets of indices specifying divisions 
  ;...of Massey's data into different populations

  if keyword_set(altrings) then begin
     div = strpos(altrings, '_', /reverse_search)
     if div NE 0 then begin
        if strpos(altrings, '_A') NE -1 then begin
           divas = strpos(altrings, '_A')
           divs = strcompress(strmid(altrings, divs+2, 10),/remove_all)
        endif
        rings = strmid(altrings, 1, div-1)
        if n_elements(divas) NE 0 then lim = divas-1-(div+2) else lim = strlen(altrings)-(div+2) 
        sky = strmid(altrings, div+2, lim) 
     endif else begin
        if strpos(altrings, 'S') NE -1 then begin
           if n_elements(divas) NE 0 then lim = divas-1-(div+2) else lim = strlen(altrings)-(div+2) 
           sky = strmid(altrings, div+2, lim) 
        endif else begin
           if strpos(altrings, 'A') NE -1 then begin
              divas = strpos(altrings, '_A')
              divs = strcompress(strmid(altrings, divs+2, 10),/remove_all)
           endif else rings = strmid(altrings, 1, div-1)
        endelse     
     endelse
     
     if n_elements(divs) NE 0 then mult = divs else mult = 1
     if n_elements(rings+sky) EQ 0 then begin
        rings = 2  & sky = 1
     endif
     nems = strarr(mult*(rings+sky))
     nems[0:mult-1] = 'core'
     if n_elements(divs) NE 0 then begin
        for c = 0, mult-1 do begin
           nems[b+c-1]+= '_'+strcompress(string(c),/remove_all)
        endfor
     endif
     
     for b = 1, n_elements(nems)/mult-1 do begin
        if n_elements(nems) GT rings then begin
           nems[b*mult:(b+1)*mult-1] = 'sky_'+strcompress(string(b-rings), /remove_all)
           if n_elements(divs) NE 0 then begin
              for c = 0, mult-1 do begin
                 nems[b+c-1]+= '_'+strcompress(string(c),/remove_all)
              endfor
           endif
        endif else begin
           nems[b*mult:(b+1)*mult-1] = 'inner_ring_'+strcompress(string(b), /remove_all)
           if n_elements(divs) NE 0 then begin
              for c = 1, mult-1 do begin
                 nems[b+c-1]+= '_'+strcompress(string(c),/remove_all)
              endfor
           endif
        endelse
     endfor     
  endif else nems = ['core', 'inner_ring', 'sky']
  if keyword_set(combo) then nems = [nems[0:1], 'combo']
                                ;getting the population names
                                ;sorted out (will deal with
                                ;combo's inflexibility later)

;want to add in something to propagate out the errors, which are
;included in members and massey as e[mag or color]

  option_label = bit 
  restore, 'massey2007.sav'  ;structure of color/mags called stars
  for i = 0, n_elements(nems)-1 do begin 
     if keyword_set(combo) AND i EQ n_elements(nems)-1 then begin
        ind = [[0,0]]
        for j = 0, n_elements(nems)-2 do ind = [[ind],[pops.(j+1)]] 
        ind = ind[*,1:*]
        inds2 = matchem(pops.key, ind)
     endif else inds2 = matchem(pops.key, pops.(i+1))
     ind = pops.(i+1)
     mv = stars.v  & mv = mv[inds2]
     mdbv = stars.dbv  & mdbv = mdbv[inds2]
     mdvr = stars.dvr  & mdvr = mdvr[inds2]
     ;obtaining the stars for the desired population

     restore, 'members.sav'
     locs = 0
     for j = 0,n_elements(members.inds)-1 do begin
        temp = where(inds2 EQ (members.inds)[j])
        if total(temp) NE -1 then locs = [locs,temp]
     endfor
     locs = locs[1:*]

     contn = fltarr(n_elements(mv))
     if keyword_set(combo) AND i EQ n_elements(nems)-1 then begin 
        fg = reform(foreground.(0))
        for j = 1, n_elements(nems)-2 do fg = [fg,n_elements((pops.(j))[0,*])+reform(foreground.(j))]
     endif else fg = reform(foreground.(i))
     fg = matchem(pops.key, ind[*, fg])
     contn[fg] = reform(info.cont_weight)
     contn[where(contn EQ 0)] = 1.
     contn[locs] = reform(info.mem_weight)
     ;marking foreground contamination and upweighting confirmed members
     
     if keyword_set(masscut) then begin
        iso_all, mv, mdbv, nems[i]+option_label, weight = contn, /masscut
     endif else begin
        iso_all, mv, mdbv, nems[i]+option_label, weight = contn
     endelse
  endfor

end


  



function run_prep, pop
;+
;PURPOSE:
;This function will create a unique tag for a given run, and generate a
;directory in which to store the results from that run

;INPUT:(None)

;OUTPUT:
;run = the unique run tag
;-
  
  run_str = systime()
  ;finding the date for the run

  yr_pos = strpos(run_str, ' ', /reverse_search) & year = strmid(run_str,yr_pos+3,2)
  mon_pos = strpos(run_str, ' ') & month = strmid(run_str,mon_pos+1,3)
  months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
  month = strcompress(string(where(strpos(months,month) NE -1)+1),/remove_all)
  if strlen(month) LT 2 then month = '0'+month
  day_pos = mon_pos+5 & day = strcompress(strmid(run_str,day_pos,2),/remove_all)
  if strlen(day) LT 2 then day = '0'+day
  run= year+month+day
  ;obtaining a run identifier on the basis of the date

  check = file_search('results/'+run+'*/*')

  if (check)[0] EQ '' then run += '_0' else begin
     run_temp=run
     run_temp += strcompress('_'+string( $
            fix(strmid(check[-1],strpos(check[-1],'/',/reverse_search)-1,1)$
                  )),/remove_all)
     if (file_search('results/'+run_temp+'/'+pop+'fitdata.txt'))[0] NE '' then begin
        run += strcompress('_'+string( $
               fix(strmid(check[-1],strpos(check[-1],'/',/reverse_search)-1,1)$
                  )+1),/remove_all)
     endif else begin
        run = run_temp
     endelse
  endelse

  ;checking to see if a run has already been performed that day
  ; and creating a later run if so
                        
  if (file_search('results/*'))[0] EQ '' then spawn, 'mkdir results'
  if (file_search('results/'+run+'/*'))[0] EQ '' OR (file_search('results/'+run+'/'+pop+'fitdata.txt'))[0] NE '' then begin 
     spawn, 'mkdir '+run
     spawn, 'mv '+run+' results'
     spawn, 'cp results/test.txt results/'+run
  ;creating a directory to save the results of a run
  endif

  return, run
end




pro iso_all, v, dbv, pop, weight = weight, masscut = masscut
;+
;PURPOSE:
;This procedure will fit every isochrone to a given population of stars
;from the Massey data

;INPUT:
;b,dbv = the mag/color of the stars that form the population
;pop = the population given as a string for the filename
;\weight = an array specifying the weighting of each star on the basis
;of the foreground contamination screening
;\masscut = by default the entire isochrone (up to the first turn) is
;used for the fit, this option uses only supergiant appropriate masses
;along the isochrone

;OUTPUT:
;-
;this function deals with the isochrones

  if keyword_set(masscut) then pop += '_mc'
  run = run_prep(pop)
  
  grid_mkr, v, dbv, offset=pop, /save, weight = weight

  openw, 110, 'results/'+run+'/'+pop+'_fitinfo.txt' 
  openw, 111, 'results/'+run+'/'+pop+'_fitdata.txt' 
  tab = ','
  printf, 111, ";  Z", tab, "                log(Age)", tab, "       E(B-V)",$
          tab, "         A_v", tab, "        Score"   
  ;keeping a log of the fit info
  
  psopen, 'results/'+run+'/'+pop+'_fits.ps', /inches, /color, /psfonts, xsize = 12, ysize = 8
  !p.font = 1 & !p.multi = [0,3,2] & loadct, 39 ;setting plot stuff
  restore, '~/ic10/Geneva/parsed/general_info.sav'
  
  isochrones = file_search('~/ic10/Geneva/parsed/*.sav')
  isochrones = isochrones[where(strpos(isochrones, 'general') EQ -1)]
  if keyword_set(masscut) then begin
     isochrones = isochrones[where(strpos(isochrones,'_mc_') NE -1)]
  endif else begin
     isochrones = isochrones[where(strpos(isochrones,'_mc_') EQ -1)]
  endelse


  for i = 0, n_elements(isochrones)-1 do begin
     pos1 = strpos(isochrones[i],'/',/reverse_search)+1
     pos2 = strpos(isochrones[i],'_')
     met = strmid(isochrones[i], pos1,pos2-pos1)
     model = strmid(isochrones[i], pos2+1,1)

     restore, isochrones[i] ;restores structure in_struct
     for j = 0, n_tags(in_struct)-1 do begin ;iterating through ages
        age = uages[j] ;getting the age for the plot and log

        iso_v = reform((in_struct.(j))[1,*])
        iso_dbv = reform((in_struct.(j))[2,*])
        iso_v = iso_v[where(iso_v NE -99)]
        iso_dbv = iso_dbv[where(iso_dbv NE -99)]
        if n_elements(iso_v) LT 10 OR n_elements(iso_dbv) LT 10 then goto, nope
        ;selecting the V magnitude and B-V color values for the isochrone

        niso = repop(iso_v,iso_dbv);,/check)
        iso_dbv = niso[0,*] & iso_v = niso[1,*]
        if total(where(iso_v EQ iso_v[1])) NE 1 then begin
           iso_dbv = iso_dbv[0:(where(iso_v EQ iso_v[1]))[1]-1]
           iso_v = iso_v[0:(where(iso_v EQ iso_v[1]))[1]-1]
        endif
        turn = loop_check(iso_v,iso_dbv);,/check)
        if turn[0] LT .75*n_elements(iso_v) AND n_elements(turn) GT 1 then turn[0] = turn[1]
;        iso_dbv1 = iso_dbv[0:turn[0]] & iso_v1 = iso_v[0:turn[0]]
        iso_dbv1 = iso_dbv & iso_v1 = iso_v;[0:turn[0]]
        ;repopulate and get up to first bend

;        if n_elements(turn) GT 1 then begin
;           !p.multi = 0  &  loadct, 39
;           plot, iso_dbv, iso_v, /xstyle,/ystyle
;           oplot, iso_dbv[turn], iso_v[turn], psym = 2, color = 75
;           stop
;        endif
;        goto, jump

        fit = isofitter(iso_v1, iso_dbv1, offset = pop)
        ;isofitter returns [hor shift, vert shift, score]
        loadct, 39

        plot, dbv, v, /nodata, xr=[-1,3], yr=[25.5,13], xtitle = "B-V", $
             ytitle = "V", title = "("+strcompress(string(age))+$
              " log(yr), Z=."+string(met)+") Isochrone Fit", charsize =1.5
        oplot, dbv, v, psym = 3
        oplot, iso_dbv[0:turn[0]]+fit[0],iso_v[0:turn[0]]+fit[1], color = 220
        oplot, iso_dbv[turn[0]:*]+fit[0],iso_v[turn[0]:*]+fit[1], color = 190
        ;plotting the isochrone to emphasize the portion that is fit
        print, fit

        ;xyouts, 'Score: '+strmid(strcompress(string(fit[2]),/remove_all),0,5)
        ;adding the score for visual comparison


        printf, 110, "File: "+ isochrones[i]

        printf, 110, "Metallicity: ", met,",  Age: ",age,",  Model: ", model
        printf, 110, "Horizontal Shift: ", fit[0], ",  Vertical Shift: ", fit[1]  
        printf, 110, "Score: ", fit[2]
        printf, 110, ""
        printf, 110, ""
        ;adding the fitting information to the human-readable log

        printf, 111, float('.'+met), tab, age, tab, fit[0], tab, fit[1], tab, fit[2]   
        ;adding info to the IDL-readable file
        jump: print, ""
        nope: if age GE 8.5 then break ;fits past this age are pointless
     endfor
  endfor
  !p.font = 0 & !p.multi = 0 & loadct, 0 ;returning to normal 
  close, 110 & close, 111
  psclose
  ;spawn, 'rm '+pop+'_massey_grid.sav'

  
end




function isofitter, iv, idbv, offset = offset
;+
;PURPOSE:
;This function will take the V magnitude and B-V color values for the
;isochrone, and the density grid defined by the data points, and shift
;the isochrone through all allowed values of shifts based on the Massey
;paper values to find the best fit.

;INPUT:
;iv, idbv = the V magnitude and B-V color values of points along the
;isochrone
;grid = the density grid from grid_mkr

;OUTPUT:(None)
;-

;need to restore the grid...actually we don't, I made grid_mkr
;awesome apparently. Woo!
;--->do actually need to make plots of the grids and the variety of values

  dist_mod = 5*alog10(660e3)-5
  ebv = findgen(20)/19.*.04-.02 + .81
  r_v = findgen(20)*1.9/19.+3.1
  ;setting up the values to shift the isochrone by

  score = fltarr(n_elements(ebv), n_elements(r_v))
  for i = 0, n_elements(ebv)-1 do begin
     for j = 0, n_elements(r_v)-1 do begin
        new_iv = iv+dist_mod+r_v[j]*ebv[i] & new_idbv = idbv+ebv[i]
        fit = 1
;SHOULD EVENIT ACTUALLY BE SET?
        if keyword_set(offset) then begin
           grid_mkr, new_iv, new_idbv, match = fit, /evenit, offset =offset
        endif else begin
           grid_mkr, new_iv, new_idbv, match = fit, /evenit
        endelse           
        temp = fit.grid
        score[i,j] = total(temp[fit.ind])
     endfor
  endfor
  ;iterating through all of the possible shifts
  
  dummy = max(score, loc)
  loc = array_indices(score, loc)
  shifts = [ebv[loc[0]], r_v[loc[1]]*ebv[loc[0]]+dist_mod, dummy]
  ;identifying what the optimal fit of the isochrone was

  return, shifts
end 



function repop, iv, ibv, check = check
;+
;PURPOSE:
;This function will take an isochrone and distribute the points along it
;evenly to ensure an unbiased fit

;INPUT:
;iv, ibv = the V magnitude and B-V color of each point on the isochrone
;/check = this keyword is to allow me to check that this function is working

;OUTPUT:
;niso = a 2 by n_elements(V) array containing the repopulated isochrone,
;first column = B-V, second column = V
;-

  niv = iv[0] & nibv = ibv[0]
  for i = 1, n_elements(iv)-1 do begin
     niv = [niv, iv[i-1]+findgen(20)*(iv[i]-iv[i-1])/19.]
     nibv = [nibv, ibv[i-1]+findgen(20)*(ibv[i]-ibv[i-1])/19.]
  endfor
  ;generating overpopulated versions of the isochrone values

  dists = sqrt((niv[1:*]-(shift(niv,1))[1:*])^2+(nibv[1:*]-(shift(nibv,1))[1:*])^2)
  dist_tot = total(dists) & interval = dist_tot/n_elements(iv)
  ;finding the distances along the isochrone so that points can be interspersed more evenly

  nv = niv[0] & nbv = nibv[0]
  ind = 0
  for i = 1, n_elements(iv)-1 do begin
     j = 1 & sum = 0 
     while sum LT interval do begin
        if j+ind GE n_elements(dists) then break
        sum += dists[j+ind]
        j += 1
     endwhile
     if j+ind GE n_elements(dists) then ind = -1 else ind += j-1
     nv = [nv, niv[ind]]
     nbv = [nbv, nibv[ind]]
     if j+ind GE n_elements(dists) then break
  endfor
  ;taking the first point that results in a relatively even distance between points

  if keyword_set(check) then begin
     !p.multi = [0,2,1]
     plot, ibv, iv, psym = 3, /xstyle, /ystyle, xtitle = "B-V", ytitle = "V", title = "Initial Isochrone"
     plot, nbv, nv, psym = 3, /xstyle, /ystyle, xtitle = "B-V", ytitle = "V", title = "Repopulated Isochrone"
     !p.multi = 0
     ;plotting to check the changes made
  endif

  niso = transpose([[nbv],[nv]])
  return, niso
end




function loop_check, iv, ibv, check = check
;+
;PURPOSE:
;This function will identify the points at which the isochrone possesses
;a loop on the basis of the changes in slope of the curve so that they
;can be dealt with in such a way that the fit is not biased by them

;INPUT:
;iv,ibv = the V magnitude and B-V values of the isochrone after repopulation

;OUTPUT:
;turn_pts = a (number of turns) array possessing the index of the ibv
;and iv arrays where a turn occurs
;-
  far = 7 ;marking the distance a "far" spot on the isochrone is a way (checking
          ;how consistent the slope of the isochrone is prior to the detected "turn")

  ;dx = ibv[1:*]-(shift(ibv,1))[1:*]
  ;xdev = where((dx/abs(dx))[1:*] NE (shift(dx,1)/shift(abs(dx),1))[1:*] AND $
  ;             (dx/abs(dx))[1:*] NE (shift(dx,far)/shift(abs(dx),far))[1:*])

  dy = iv[1:*]-(shift(iv,1))[1:*]
  ydev = where((dy/abs(dy))[1:*] NE (shift(dy,1)/shift(abs(dy),1))[1:*] AND $
               (dy/abs(dy))[1:*] NE (shift(dy,far)/shift(abs(dy),far))[1:*])
  
  spots = ydev ;uniq([xdev,ydev], sort([xdev,ydev]))+1

  if keyword_set(check) then begin
     !p.multi = 0
     loadct, 39
     plot, ibv, iv, /xstyle, /ystyle, xtitle = "B-V", ytitle = "V", title = "Isochrone"
     oplot, [-99,ibv[spots]], [-99,iv[spots]], psym = 2, color = 250
     loadct, 0
  endif

  return, spots
end




pro grid_mkr, v, dbv, match = match, save = save, evenit = evenit, weight = weight, offset = offset
;+
;PURPOSE:
;This procedure will generate a grid with the density of stars in a
;particular region of color-magnitude space, and identify which stars
;belong to each grid point

;INPUT:
;v, dbv = the V magnitude and B-V color
;/match = This keyword will take a grid defined on one data set (say the
;Massey data), and find the location of data points from another data
;set (say our data), upon it
;/save = This keyword will save the grid as massey_grid by default
;/evenit = This keyword removes points in the grid below 1 sigma above
;the average value of the grid and sets all non-zero points to the same value
;/weight = optional input to allow not all stars to be counted as equal
;in designing the grid
;/offset = an optional input that prevents restoring issues when running
;the fit in parallel by creating a unique massey_grid file

;OUTPUT: 
;-

  if keyword_set(match) then begin
     if keyword_set(offset) then begin
        restore, offset+'_massey_grid.sav'
     endif else restore, 'massey_grid.sav'
     
     if keyword_set(evenit) then begin
        grid[where(grid LT median(grid)+stddev(grid))] = 0
        grid[where(grid NE 0)] = min(grid[where(grid NE 0)])
     endif

     locs = v*0
     for i = 0, n_elements(var)-2 do begin
        for j = 0, n_elements(dbvar)-2 do begin
           inds = where(v LT var[i+1] AND v GE var[i] AND $
                        dbv LT dbvar[j+1] AND dbv GE dbvar[j])
           if total(inds) NE -1 then locs[inds] = i*(n_elements(dbvar)-1)+j
        endfor
     endfor

     ;finding the positions of data points on the Massey grid

     grid_info = {grid: grid, ind: locs}
     match = grid_info
     goto, done
  endif

  if keyword_set(weight) then weight = weight else weight = fltarr(n_elements(v))+1.

  vr = (max(v)-min(v)) & dbvr = (max(dbv)-min(dbv[where(dbv GT -1)]))
  if vr LT dbvr then begin
     var = findgen(25*vr/dbvr)*vr/(25*vr/dbvr-1.)+min(v)
     dbvar = findgen(25)*dbvr/24.+min(dbv[where(dbv GT -1)])
     grid = fltarr(25-1,25*vr/dbvr-1)
  endif else begin
     var = findgen(25)*vr/24.+min(v)
     dbvar = findgen(25*dbvr/vr)*dbvr/(25*dbvr/vr-1.)+min(dbv[where(dbv GT -1)])
     grid = fltarr(25*dbvr/vr-1,25-1)
  endelse
  ;creating the ranges and arrays for the grid

  grid_locs = v*0
  for i = 0, n_elements(var)-2 do begin
     for j = 0, n_elements(dbvar)-2 do begin
        inds = where(v LT var[i+1] AND v GE var[i] AND $
                     dbv LT dbvar[j+1] AND dbv GE dbvar[j])
        grid[j,i] = total((inds NE -1)*weight[inds])
        if grid[j,i] NE 0 then grid_locs[inds] = i*(n_elements(dbvar)-1)+j
     endfor
  endfor
  ;creating the grid and finding the positions in the grid where each data point falls

     
  if keyword_set(evenit) then begin
     grid[where(grid LT median(grid)+stddev(grid))] = 0
     grid[where(grid NE 0)] = min(grid[where(grid NE 0)])
  endif

  if keyword_set(save) then begin
     if keyword_set(offset) then begin
        save, var, dbvar, grid, grid_locs, filename = offset+'_massey_grid.sav'
     endif else begin
        save, var, dbvar, grid, grid_locs, filename = 'massey_grid.sav'
     endelse
  endif
  done: if 3 EQ 4 then print, "huh?"
end




