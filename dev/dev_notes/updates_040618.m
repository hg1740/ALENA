% Tweaks / fixes
% 
%
% XML Export
%
%  a) Added ability to export tables to xml in neater fashion (one less level of unnecessary hierarchy that
%      had crept in as part of my first attempt).
%
%      Tables exported by recursively calling the export routine on each column of the table in turn.
%
%  b) Added ability to export logical and datetime types, which I seem to have forgotten about previously.
%
%      This tweak highlighted a particular issue with the "Locked" flags that control user's ability to
%       modify structure or parameters associated with a model.  If a file is exported with any locked flag
%       raised, it renders the resultant file unimportable, because the import process relies on being able
%       to reconstruct the model piece by piece as the XML is parsed.
%
%      Addressed this issue by allowing an option to be passed in to the "add" method used to build the model,
%       allowing state of the Locked flags to be ignored, for the duration of that particular "add" operation.
%       The import process then can be treated as special, by specifying this new "ignoreLockedFlags" option
%       each time it makes a call to "add".  A bit of a pain, but it works ok.
%
%      An alternative approach would be to exclude these properties from the export file,
%       see item (e) below - but not sure if we'd want this or not for these properties - to discuss.
%
%  c) Added ability to export cell-arrays to xml properly (was previously limited to just cellstr, which are easier to deal with).
%
%      Cell export now implemented in similar fashion to table export, with a recursive call to each element in turn.
%      One difference however is that you have to think of a name using which to tag the content,
%      unlike the table which has a natural column name you can use.  Would be noce to use e.g. {idx},
%      where idx is the index into the cell, so the XML would look "readable" to someone familiar with cell array syntax,
%      but {} characters appear to be illegal in xml tags, so I've just labelled each value "_idx" instead.
%
%  d) Non-scalar quantities exported to xml now include their "size" as an attribute,
%      (a more generic approach than exporting number of rows and/or columns explicitly).
%      should make it much easier to reconstruct content correctly during import.
%
%  e) Exporting the "BeamModel" class raises issues also previously unhandled in the mvc baseclass,
%      including export of structure arrays, and of n-d arrays where n > 2.
%
%      These now handled, although I'm not sure as well as they could be (again by using extra hierarchy
%      with the "_idx" tag used to identify element idx of a structure array, or page idx of an n-d array)
%
%  f) This additional capability of the export routine now means that an AWI session is exported in fuller
%   detail than previously, now including content such as the Audit Trail, which I'm not sure is really required
%   in an export - so maybe we now need a means of down-selecting what is actually exported.
%
%  One option might be to link it to the "Nameable" class, which manages the different tab pages in the
%   Properties view - this already allows for attributes that control whether a property is visible,
%   or enabled for editing.  I've added an attribute to each property group, or each property, controlling
%   whether or not it is exported, and adjusted the constructor for awi.model.Component, to configure it
%   so marker properties are only exported if the "Export" attribute is true.  The attribute defaults to true,
%   so it is included in export unless you say otherwise.
%
% XML Import
%
%   a) Seems to me like the "Reimport..." menu item would be more useful if it told you the
%       name of the file that it would use, if invoked.  That then sounds similar to the
%       functionality of the "File->Recent" menu, with multiple files on a "most recently used" list.
%
%      Have therefore introduced the concept of a "most recently imported" list,
%       therefore consistency across "open" and "import" behaviour (and shared code),
%       whilst providing for a more user friendly interface at the same time.
%
%   b) The import progress dialog is presented as we've seen previously, but when using an
%       extended desktop it appears on the default monitor, not necessarily the one on which
%       the framework GUI is displayed.
%
%      Have tweaked it so it is initialised centred on the GUI from which it was invoked.
%
%      Have taken the opportunity to make this dialog a generic method, implemented in
%       mvc.mixin.UiTools, rather than anything too specific to the Import process, as I think
%       it will be useful elsewhere too (e.g. for displaying progress during analyses that iterate,
%       or seek to minimise something till they converge.
%
%      Something to consider:  should the import progress window actually be treated as another "view"
%       into the application object ?  At present, it is displayed as a separate, modal, figure window.
%       Could be re-implemented as an "Import Progress View", or it could just simply be the Audit Trail view ?
%
%  c)  Previously, if the import process fails, the progress dialog was closed and then the
%       failure message was displayed in a separate error dialog.
%       
%      Have tweaked this, so that any failure messages are written to the progress dialog,
%       which remains open until the process is complete.
%
%      This is much more user-friendly, for example the BUG file now imports the aircraft,
%       and reports that the "loadcases" file referenced within the BUG file could not be found.
%       At least you get to see what was importable as far as the error.
%       (and this particular problem with importing files that reference other files has been
%        addressed - see point (d)).
%
%      This also makes it easy them to retain the list of import messages, and store it in the Audit Trail.
%
%  d) The import file folder name is now passed down through the recursive calls of the import routine,
%      so one file imported by reference from another file should now work, if both files are in same folder.
%
%  e) Import routine now aware of additional data types supported by export (logical, datetime, cellstr),
%      as well as structure arrays and n-d arrays.
%
%  Other Export / Import issues
%
%  a) The generic (i.e. mvc) exportable/importable classes support XLS format as well as XML.
%     Have tweaked the XLS export to make it slightly better, and XLS import to make use of the
%     progress window in similar fashion to XML.
%
%     However, there's a bit of an issue with how best to represent data structures, tables
%     and n-d arrays in an XLS export - nothing that could not be handled, but it would take a bit
%     of effort, and not strictly required for AWI (I guess) so perhaps best left for time being.
%
%  b) TODO: change to the import workflow also needs to also be reflected in the import_fm4 method,
%      which I can't test at my end, so something to work on at next session at UoB.
%
%  c) Exporting a file in mat-file format is supported, but the implementation basically passes the call
%     on to "save", so the mat-file has the same content and structure as if the session had been
%     saved (i.e. the exported mat-file is basically the same as a saved awi-file).
%
%     This has unexpected consequence - the mat-file then appears on the "recently saved" list,
%     because the code that maintains this list doesn't know whether it was arrived at by means of save,
%     or export.  Have provided a "-nomru" option to save/open methods of Serializable, so this can be suppressed.
%
% Results Viewer
%
%   a) I figured out a way of dispaying the legend in the vertical toolstrip on RHS of the view, rather
%      than within each chart themselves (the way I have done it is completely clunky - I had to clone
%      one of the charts, parent the clone in the vertical toolstrip, add a legend to it, then hide it,
%      very unelegant, but bottom line is that it seems to work).
%
%   b) Easy to add a custom hit handler to each item displayed in a legend, so I have put in place a
%      callback that toggles the visible state of the lines in all charts in the view, when you click
%      on the corresponding legend string.  Rather neat, I think.
%
%      TODO: Still need to figure out how to achieve same thing programatically, on basis that there should
%      be nothing you can do through the UI that you can't do from a script.
%
%   c) Have added a rudimentary csv export facility.
%
%      TODO: customise for Airbus-specific format / file structure, if necessary.
%
%
% Layered Drawing view – I noticed a couple of features at our last meeting that were not quite what I had intended 
%
%   a)	Selection list offers you items that are not Drawable, which is not really appropriate – needs down-selection
%        before the popup menu is populated.
%
%       It turns out that the list IS down-selected to Drawable items.  The issue is that things like LoadCase,
%        which you'd not expect to be drawable, are derived from Entity, which is in turn a DrawableThing.
%
%       TO CONSIDER: Is this what we actually want ?  I think I recall Chris S had a good reason for making them Drawable....
%        But for tinme being, leaving them as Drawable but setting Visible = false achieves the desired fix
%        to the LayeredDrawing selection list.
%   
%   b)	Adding a new object to the collection that is “Drawable” should simply “just work”,
%        i.e. it should be offered immediately as an item in the drawing view, and as long as “drawElement” 
%        is implemented it should draw itself.
%
%       It seemed to very nearly work at last meeting but not quite – need to trace through the code more 
%       carefully (i.e. without you lot watching) to better understand why it didn’t quite work as intended.
%
%       Having done so - I can't find the sequence of events that makes it break, need to look into further
%        at next session.
