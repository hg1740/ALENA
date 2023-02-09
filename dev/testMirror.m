function c = testMirror
%testMirror Test the 'awi.mixin.Mirrorable' functionality.

close all force

c(1) = awi.model.Component('Name', 'Component(1)');
c(1).XOffset = 1;
c(1).YOffset = 2;
c(1).ZOffset =  3;

c(2) = awi.model.Component('Name', 'Component(1) (Mirror)', '-Mirror', c(1));

draw(c);

end

