
% TODO - Round up number instead of truncation

function OutStr = num2_8ch(InpNum)

Field  = 8;
Field2 = 8;

InpStr = num2str(InpNum,Field-1);
LenInp = length(InpStr);

I  = find(InpStr == '.');
Ie = find(InpStr == 'e');

if isempty(I) && isempty(Ie)
    
    if LenInp >= Field
        OutStr = num2str(InpNum,Field-2);
        LenOut = length(OutStr);
        
        I  = find(OutStr == '.');
        Ie = find(OutStr == 'e');
        
        Res = LenOut - Field;
        OutStr(Ie-Res:Ie-1) = [];
        
    else
        OutStr = [InpStr '.'];
    end
    
elseif ~isempty(I) && isempty(Ie)
    
    if LenInp > Field
        InpStr(Field + 1:LenInp) = []; % may need to round this instead
    end
    
    OutStr = InpStr;
else
    if LenInp > Field
        Res = LenInp - Field;
        InpStr(Ie-Res:Ie-1) = [];
    end
    
    OutStr = InpStr;
end

end
