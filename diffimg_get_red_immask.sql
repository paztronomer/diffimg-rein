-- Query to get new arrived and processed data
-- The below returns filenames, path 
-- Both works!

select distinct(fai.filename), fai.path, fai.compression
    from file_archive_info fai, desfile d, pfw_attempt_val val
    where d.pfw_attempt_id=1153411
    and val.key='expnum'
    and to_number(val.val,'999999')=615518
    and d.filename=fai.filename
    and d.filetype='red_immask'
    order by fai.filename;

select fai.path, fai.filename, fai.compression
    from image im, file_archive_info fai
    where im.pfw_attempt_id=1153411 
    and im.filetype='red_immask' 
    and im.expnum=615518 
    and fai.filename=im.filename;
