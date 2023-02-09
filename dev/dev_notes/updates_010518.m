% Tweaks / fixes
% 
% 1.  Accelerators added to File>New etc, so you can now use "Ctrl-N",
%     "Ctrl-O", "Ctrl-S", "Ctrl-W" as shortcuts
%
% 2.  The "File->Close" function was broken in last version (it did nothing,
%     because of an inadvertent "return" statement I had put in for some reason
%     I now forget.
%
% 3.  Added an option to "File->Close All", so you can more easily close multiple
%     sessions with a single command from the GUI
%
% 4.  Size of "HotSpot" in View Manager - AT LAST, cracked it.  View Manager
%     context menu now available on right-click on the actual tab, rather than
%     in the region between the tab and the content, which can become tiny when
%     the view is fairly busy.  It only took 2 lines of code.
%
% 5.  Analyse menu: now honours the tokens '>' indicating sub-menu, '|'
%     indicating a separator (and '!' indicating an accelerator), 
%
% 6.  You might have noticed (or quite possibly not) that on invoking "Edit>Duplicate>This"
%     on either the top-level Framework node, or the Aircraft node, results in an error
%     because, quite correctly, the tool is telling you that there is no valid place within
%     the framework session that a duplicate of either item could be stored.
%
%     But it seems to me to be entirely reasonable that the tool should allow these items to
%     be duplicated - but doing so results in the duplicate being sent to a new instance of
%     the framework, rather than the current one.  Behaviour tweaked accordingly.
%
%     Furthermore, if you invoke "duplicate this" on a load case or result set, you get an
%     unexpected choice to two destinations for the duplicate, because of some slightly dodgy
%     logic in my generic code that works out where the duplicate could be parented.  Now fixed.
%
%     Also, I'd let a bug creep into "Edit->Duplicate->Children...", suspect it hasn't
%     worked for a while.  Now fixed.
%
% 7.  Found a bug in mvc.mixin.Serializable, which would (rather bizarrely) cause an error if you try to
%     load / save session whilst "rotate" is switched on in the axes (which means it doesn't get picked
%     up by my auto test script).  Now fixed.
%
% 8.  Tree node expand / collapse state now retained during interactive edit.
%
% 9.  Bug fix in mvc.mixin.Searchable
%
% 10. Rather crude, but it proved straight-forward to add a feature whereby the aircraft
%     "Trim" can be set by picking from a list of available TrimResults.
%     The RotMat associated with the selected TrimResult is set as the "Orientation" of the aircraft,
%     so it redraws itself as specified by the trim analysis.
%
%     I appear to have the aircraft pitching down, when I'd expect it to pitch up, so there's
%     probably a sign issue somewhere, which we can sort out if you think this feature is worthwhile.
%     The next logical step would be to set the control surface deflections in the model
%     from those picked up from the trim result - but I'm not sure we're quite in a
%     position to do so yet, I'm not sure how to interpret the 'ControlDeflection' field
%     in the TrimResult object.
%
% 11. Framework was far too easily broken simply by not having a Cruise loadcase - because we
%     had set it up such that the very existence of the SizeAnalysis view was deemed invalid
%     in the ansence of a cruide loadcase.  This is far too onerous - better instead to allow the
%     view to be created, and only report the missing loadcase when the "analyse" button is clicked.
%
%     I then got fed up importing the sugar_volt aircraft, running "size" to see it error immediately
%     because it alleged there was no Cruise load case, when in fact there is a manouevre case with
%     loadfactor == 1, so I made the analysis tolerable to this combination of inputs.
%     Easily reverted if you think that is a step too far.



