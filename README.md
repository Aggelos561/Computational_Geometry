                    -- Εργασια 1 --

-- Πολυγωνοποίηση σημειοσυνόλου με την χρήση της βιβλιοθήκης CGAL --

Άγγελος Δωρόθεος Χατζόπουλος 1115201900217

Γρηγόριος Μουλκιώτης 1115201900117


-- polygonization.hpp - polygonization.cpp --
Αρχικά ορίζουμε την κλάση Polygonization που περιέχει τον τρόπο που επιλέγονται οι ακμές στην επιλεγμένη μέθοδο(edgeSelection),το εμβαδόν και το ratio(totalArea,ratio αντίστοιχα) περιέχει επιπλέον μεθοδους για επέκταση της πολυγωνικής γραμμής(expandPolygonLine),την διαγραφή ευθύγραμμων τμημάτων(deleteSegment),την εισαγωγή σημείων στο πολύγωνο όταν αυτα βρίσκονται στο κυρτό περίβλημα(forceInsertPoint),την διαγραφή ευθυγράμμων τμημάτων(deleteSegment) και getters και τις μεταβλητές της κλάσης. 


-- incremental.hpp - incremental.cpp --
Για τον αυξητικό αλγόριθμο χρησιμοποιούμε την κλάση Incremental που κληρωνομεί απο την κλάση Polygonization η οποία παρέχει μεθόδους για αρχικοποίηση του τριγώνου(initializeTriangle),για την εύρεση των κόκκινων,ορατών και επιλογή ορατής ακμής(chooseVisibleSegment,getRedSegments,findVisibleSegments),διάταξη σημείων με βάση τις συντεταγμένες τους(sortYAsc,sortYDesc,sortPoints) όταν μας δοθεί ως όρισμα απο το command line ο τρόπος αρχικοποίησης του κυρτού πολυγώνου και αποθηκεύεται στην μεταβλητή της κλάσης initialization. Η μέθοδος start στην κλάση incremental καλείται απο την main για την εκτέλεση του αλγορίθμου.


-- convexHull.hpp - convexHull.cpp --
Για τον αλγόριθμο κυρτού περιβλήματος χρησιμοποιούμε την κλάση convexHull που κληρωνομεί απο την κλάση Polygonization και παρέχει μεθόδους για την αρχικοποίηση του κυρτού περιβλήματος(initializeConvexHull),την εύρεση του κοντινότερου σημείου απο κάθε ακμή και επιλογή του καταληλλότερου για εισαγωγή στο πολύγωνο(findBestPoint,findVisiblePoints,insertBestPoint) και υπολογισμό του εμβαδού(calcArea). Η μέθοδος start στην κλάση convexHull καλείται απο την main για την εκτέλεση του αλγoρίθμου.


-- dataio.hpp - dataio.cpp --
Για το io χρησιμοποιούμε ενα namespace που περιέχει συνάρτηση για διάβασμα των command line arguments(getParameters),το διάβασμα των σημείων απο αρχείο(readPoints) και την εξοδο σε output file(createResultsFile).


-- main.cpp --
Οπότε στην main απλως επιλέγουμε την κατάλληλη κλάση ανάλογα τον αλγόριθμο και επιλέγουμε την μέθοδο start για το τρέξιμο του αλγορίθμου.



Για να γίνει build απαιτείται στο directory ./src/ να τρέξουν οι εντολές:

• cgal_create_CMakeLists -s to_polygon

  --> Με την εκτέλεση της παραπάνω εντολής θα πρέπει να έχει δημιουργηθεί αυτόματα ενα αρχείο CMakeList.txt που θα περιέχει μέσα τις ακόλουθες εντολές:
    - include_directories( BEFORE ../include )
    - add_executable( to_polygon  convexHull.cpp dataio.cpp incremental.cpp main.cpp polygonization.cpp )
    - add_to_cached_list( CGAL_EXECUTABLE_TARGETS to_polygon )
    - target_link_libraries(to_polygon PRIVATE CGAL::CGAL )


• cmake -DCMAKE_BUILD_TYPE=Release .


• make


Παραδείγματα εκτέλεσης προγράμματος:

    -- Για τον αυξητικό αλγόριθμο --
    
 •   ./to_polygon -i stars-0000200.instance -o results_file.txt -algorithm incremental -edge_selection 1 -initialization 1a

 •   ./to_polygon -i stars-0000200.instance -o results_file.txt -algorithm incremental -edge_selection 3 -initialization 2b


    -- Για τον αλγόριθμο με βάση το ΚΠ --

 •   ./to_polygon -i stars-0000200.instance -o results_file.txt -algorithm convex_hull -edge_selection 1

 •   ./to_polygon -i stars-0000200.instance -o results_file.txt -algorithm convex_hull -edge_selection 2


    -i το input file (.instance με κατάλληλη δομή με tabs)
    -o το output file που θα εκτυπωθούν τα δεδομένα όπως αναφέρει η εκφώνηση
    - algorithm για την εκτέλεση συγκεκριμένου αλγορίθμου (incremental or convex_hull)
    - edge_selection για την επιλογή κατάλληλης ακμής/σημείου (random ή πρόσθεση ελάσχιστου/μέγιστου εμβαδού)
    - initialization για τον αυξητικό αλγόριθμο για την επιλογή του sorting των σημείων (1a , 1b, 2a, 2b)
