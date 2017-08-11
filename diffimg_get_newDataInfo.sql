-- Query to get new arrived and processed data
-- The above only returns info from processing assessments, not path on disk

with z as (
    select fcut.expnum, max(fcut.lastchanged_time) as evaltime
    from firstcut_eval fcut
    where fcut.analyst!='SNQUALITY'
    group by fcut.expnum
    )
select e.expnum, e.nite, e.airmass, e.obstype, e.date_obs, e.mjd_obs, e.telra, 
    e.teldec, e.radeg, e.decdeg, e.band, e.exptime, 
    val.pfw_attempt_id,
    fcut.fwhm_asec, fcut.t_eff, fcut.skybrightness
    from z, exposure e, firstcut_eval fcut, pfw_attempt_val val 
    where e.obstype='object' 
    and e.nite between 20170201 and 20170203
    and fcut.expnum=z.expnum
    and fcut.expnum=e.expnum
    and fcut.lastchanged_time=z.evaltime
    and fcut.program='survey'
    and fcut.accepted ='True'
    and fcut.processed ='True'
    and val.key='expnum'
    and to_number(val.val,'999999')=e.expnum
    and val.pfw_attempt_id=fcut.pfw_attempt_id
    order by e.nite, e.band; -- drop semicolon inside python

