% Tweaks / fixes
%
% 1. Entry-point function (AWI.m): Have cleaned this up a bit, specifically:
%    i) Include the name of the entry=point function itself in the call to Framework,
%        (needed by features inherited from mvc.mixin.Deployable)
%   ii) Include a copyright statement (or any other text you want),
%        this text is pulled through to the application "About" box
%
% 2. On starting up AWI framework, the view of an "empty" session is not helpful.
%    Have tweaked behaviour so if there are no drawable children it displays a (hopefully)
%    helpful text message, instead of the the seemingly meaningless arrows.
%
% 3. File->reimport: intended behaviour is that the import method is reinvoked on whatever file was
%    last imported.  Need to cater for when there is no last imported file - options are either to
%    suppress the 'reimport' method, or to have behaviour revert to 'import'.
%
% 4. On right-click, Edit->Add... it is slightly annoying to have to pick what you want to add
%    using a separate helper dialog. Slicker to dynamically extend the "Add" submenu with whatever is addable.
%
% 5. I found it desirable under some circumstances to be able to display properties view of two
%    different items in the tree simultaneously (for example position two "Properties" views one above
%    the other).  But the existing Framework does not allow this, because the "Selected" item in a
%    Properties view is linked to that in the Tree view (the idea being that if you change one,
%    then the other changes too).  Seemed reasonable at the time.
%    But the consequence of this is that two Properties views are doomed to always display properties
%    of the same item in the tree - because changing the selection in one Properties view is propagated
%    automatically back to the Tree view, but changing the selection in the Tree is then propagated
%    automatically onwards to the other Properties view.
%    Easy enough to optionally break this link with a new checkbox "linked" (default == true) in the
%    Properties view.
%
% 6. Same comment applies to the LayeredDrawing view, which offers the user the flexibility to draw
%    only part of the aircraft.  If this view is useful, then it is also necessary to allow the
%    selection in this view to be "unlinked" from the selection in the overall framework,
%    so you can have multiple such views, each showing a different part of the aircraft.
%
%    I also modified the layered drawing view, so it gives you options not just to view the entire
%    aircraft, or the selected node in the aircraft, but a third option of the selected node plus
%    all its children - which seems to me to be jolly useful.
%
%    Also, it is straight forward to make the drawing interactive, so you can click on a graphical feature
%    and that triggers a change in the selected item in the tree, to select whichever component of the aircraft
%    owns that graphical feature.  More intuitive than having to go back to the tree to select the thing
%    you want to interact with.
%
%    With these modifications in place, I think the "Layered Drawing" is so much more useful, when viewing
%    a model that comprises a hierachy of collectable things, that it should be the default Drawing view.
%
%    I'd also like to suggest that the "tags" applied to different bits of the lifting surface view
%    (aero panel and aero profile) are not also tagged with the object name.  The intent of the tag,
%    and the way it interacts wiht the layered drawing view, is to allow you to down-select what is drawn
%    based on its type, rather than its name.
%
% 7. I'd like to propose a change to the way that "AdditionalProperties" are handled,
%    as and when they are found during XML import.  MATLAB has built in support for such a situation,
%    which in our haste to implement this feature I overlooked, so didn't exploit.
%    It is easy to build the use of "dynamicprops" into the class hierarchy, with no apparent change
%    to the end user.
%
% 8. As it stands, the only reason the XML files import at present is because we have included a
%    "relaxation" of the schema, so that any unidentified objects defined in the XML are cast to
%    generic "Entities", and added to the tree.   (The "Mesh" item is a case in point)
%
%    Is that a good thing ?  Would it not be better to enforce a stricter regime, for example:
%    i) unidentified content is assumed to invalidate the entire file, and the import fails
%   ii) unidentified content is ignored, and the import carries on regardless
%  iii) unidentified content is put to one side, whilst the import proceeds, so the user can
%       be told about it, and decide what to do at their leisure.
%
%    I haven't changed anything in the code - this comment is just food for thought.
%
% 9. We noticed at last meeting that there was quite a lot of common code shared by allegedly
%    different methods, specifically import_xml and import_fm4.  This suggests to me that
%    what we really need is an AWI-specific version of the generic "import" routine, in which
%    the common code would be contained, with the calls to specific import routines made as
%    required, minimise code duplication.
%
%    The nature of the shared code is, I think, around what to do with content after import,
%    whatever file format it comes from.  This is what I've implemented:
%
%   i) If an imported file contains an Aircraft, it is assigned to the AWI session.
%      Since the framework is configured to allow only a single aircraft per session,
%      any existing Aircraft model must be overwritten.
%      If the application "dirty" flag is set, the user is prompted before this overwrite.
%
%  ii) If an imported file contains one or more Load Cases, they can either be assigned
%      to the current AWI session, overwriting any existing load cases, or they could be added.
%      It's not obvious to me whether one or the other of these should always happen,
%      so I ahve set the framework up to ask the user.
%
% iii) An interesting further option is if the imported file contains part of an aircraft,
%      e.g. a wing.  Surely reasonable to treat it like a loadcase, and offer the oppotunity
%      to append it to part of the existing aircraft, so you could pull in something from a
%      component library, a set of wing options for example.
%
%    This logic is in the Framework "import" method, rather than the "import_xml" method, so
%    the behaviour will be the same regardless of which file format was selected during import.
%    This has the added benefit that any Framework-specific calls, like that to Chris S' "buildable",
%    can be included in the Framework import method, so appear only once in the code, and it will
%    always run after import, regardless of file type.
%
% 10. Session persistence (i.e. File->Save, Load etc) did not quite work, but only because some
%     of the classes in the awi package were not declared "ConstructOnLoad".  Easily fixed.
%     It then is a simple matter to make the persistent session file (with the "awi" file extension)
%     also treatable as "importable" in same fashion as XML or Fame4.  Thill will be handy if you want,
%     for example, to transfer just the load cases, or just the aircraft, between two session files.
%     By treating it as an import process, rather than save/load, we can pick and choose which bits of
%     one session to transfer to the other.
%
% 11. I think it is a really good idea if all objects can be incrementally added to the framework,
%     without throwing errors due to their being partially populated.
%     Specifically, if you add an "empty" lifting surface etc, taking default properties as assigned
%     during contruction, the framework needs to not fail on unhandled errors because fields are empty.
%     That sort of behaviour makes automatic scripting of tests much harder.
%     So I have tweaked some AWI classes (ControlSurface, BluffBody, LiftingSurface) so that they
%     return empties, rather than erroring, if content has not been set.
%
% 12. Along similar lines, it is a v good idea if one can save/load or export/import a session
%     programatically (good for automated testing), but that in turn means some rules have to be
%     adhered to, such as any value that you "get" from a property is a valid value to then "set"
%     on the same property (else you risk unexpected errors being thrown during automated testing).
%     Specifically, the "ActiveSet" on buildable - the default value is "none", but this is not a
%     valid value to put back using "set". I tweaked the set method to permit this.
%
% 13. Chris S' "buildable" class - various tweaks incorporated:
%
%    i) Have changed the way Chris was invoking "build", which was by means of a line of
%       AWI-specific code inserted into mvc package.  Such "code leakage" across packages
%       is not good practice - but very easily sorted, as we can create a new listener for
%       the Framework's "onModelChanged" event, which gives exactly the same behaviour but
%       without the need for "leakage" of code between packages.
%
%       But doing so is not complete plain sailing - need to worry about execution order,
%       as "onModelChanged" now triggers separate listeners to "build" and refresh all the views.
%       Unfortunately we can't guarantee the order in which these events are actioned, so there
%       is a risk that the "build" action occurs second, which is the wrong way round.
%
%       (This is one reason why standard property listeners support "PreSet" and "PostSet" events)
%
%       I think best solution is to declare two events on the Collectable class: ModelChanging,
%       and ModelChanged.  Then we can listen for "ModelChanging" and respond with the call to "build",
%       whilst ModelChanged triggers a refresh of all the views.  Seems to work nicely.
%
%   ii) With this approach in place, I can't see the need any longer to "pause" and then "resume"
%       the listeners during the import process, and xml-import now becomes completely generic,
%       with nothing AWI-specific at all as far as I can see.
%       (Pause /resume still required during the hiarerchical call to "build")
%
%  iii) I've added an attribute "BuildStatus" to the Buildable class, then wrapped the
%       "buildElement" code in a try-catch, so if anything fails the error message can
%       be captured.  This is good practice because "buildElement" retrieves and evaluates
%       function handles, over which the Buildable code has no control, so it is a good
%       idea to allow it to fail as gracefully as possible.
%
%       If buildElement does not fail, then I write meaningful message to BuildStatus,
%       so you can see what happened at the time of the last build for any item in the collection.
%
%   iv) The BuildStatus string is now displayed in the "Build Set" tab of the properties view.
%       I've also changed the "Active Set" and "Active Set Mode" selections into popups.
%       I've also configured the propertis view so that it does not show property pages
%       associated with build sets other than the Active one.
%
% 14. Did I correctly note down in my logbook that someone said, during our meeting with Airbus
%     back in 2017, that they wanted to use the AWI tool in R2015b ?
%     I found a number of constructs in the code that were not compatible - search the codebase
%     for "verLessThan('9.3')" to see where changes were necessary.  In summary:
%
%    i) Chris S had made (limited) use of "string" variables, which did not exist in 2015b.
%   ii) My mvc package makes reference to "timetables", which did not exist in 2015b.
%       (The AWI code doesn't use them anyway, but 2016b doesn't know that.)
%  iii) The "LiftingSurface" class does some element-wise maths using implicit "repmatting"
%       across rows of a matrix - which was not supported in 2015b.
%   iv) I found one very strange behaviour that I can't explain other than a bug in MATLAB's
%       OO infrastructure in 2015b that has been fixed in 2017b.
%
%     With these issues addressed, the code runs (as far as I can tell) equally well in 2015b and 2017b.
%