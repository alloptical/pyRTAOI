function ExportPeakTrigAvgToBaseWorkspace(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
ud = getappdata(thisfig,'data');

[sel,ok] = listdlg('ListString',ud.names,'PromptString','Choose data to export','Name','Export');
if ok == 0;
    return
end

for sidx = 1:length(sel);
    eval(['op.',ud.names{sel(sidx)},'=ud.data{sel(sidx)};']);
end

thisname = get(thisfig,'Name');
defname = ['STA' strrep(thisname(1:strfind(thisname,':')-1),'-','_')];
answer = inputdlg('Choose name of variable to create','Variable Name',1,{defname});
if isempty(answer);
    return
end
assignin('base',answer{1},op);
