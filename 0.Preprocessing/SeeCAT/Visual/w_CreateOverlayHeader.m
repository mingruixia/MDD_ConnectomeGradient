function OverlayHeader=w_CreateOverlayHeader(OverlayFileName, PMax, PMin, NMin, NMax, Alpha, CBarString)
[Data, Vox, Header]=y_ReadRPI(OverlayFileName);
Header=w_ExtendHeader(Header);
Header.Vox=Vox;

Header.IsSelected=1;
Header.Raw=Data;
Data = Data .* ((Data < NMin) + (Data > PMin));
if NMax >= 0
    Data(Data<0) = 0;
end
if PMax <= 0
    Data(Data>0) = 0;
end
Header.Data=Data;
Header.NMax=NMax;
Header.NMin=NMin;
Header.PMin=PMin;
Header.PMax=PMax;
Header.cbarstring=CBarString;
Header.numTP=1;
Header.curTP=1;
Header.Transparency=Alpha;

if ischar(CBarString)
    if CBarString(end)=='+' || CBarString(end)=='-'
        PN_Flag=CBarString(end);
        CBarString=CBarString(1:end-1);
    else
        PN_Flag=[];
    end
    cbar=str2double(CBarString);
    if isnan(cbar)
        ColorMap = colormap(cbarstring);
    else
        ColorMap = y_AFNI_ColorMap(cbar);
    end
    Header.ColorMap=ColorMap;
else
    ColorMap=CBarString;
end

ColorMap = y_AdjustColorMap(ColorMap,...
    [0.75 0.75 0.75],...
    NMax,...
    NMin,...
    PMin,...
    PMax,...
    PN_Flag);
Header.ColorMap=ColorMap;

OverlayHeader=Header;