
function detrend_out=detrend_nan(indata)

    teme_data=indata;
    ind_useful=find(~isnan(teme_data));
    teme_data(isnan(teme_data))=[];
    teme_data2 = detrend(teme_data);
    detrend_out=nan(length(indata),1);
    detrend_out(ind_useful)=teme_data2;

end