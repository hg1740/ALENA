function offset = get_midline_offset(NODE1, NODE2, NODE3, OFFSET1, OFFSET3)

x1 = NODE1 + OFFSET1;
x3 = NODE3 + OFFSET3;

midnode = mean([x1 ; x3], 1);
% offset2
offset = midnode - NODE2;

end