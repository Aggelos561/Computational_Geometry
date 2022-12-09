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


• Εκτέλεση του προγράμματος όπως παρακάτω


Παραδείγματα εκτέλεσης προγράμματος:


    -- Για τον αυξητικό αλγόριθμο --
    
 •   ./to_polygon -i euro-night-0000050.instance -o results_file.txt -algorithm incremental -edge_selection 1 -initialization 1a
     Algorithm: incremental_edge_selection_1_initialization_1a
     Area: 45615704
     ratio: 0.697328
     Construction time: 0

 •   ./to_polygon -i london-0000100.instance -o results_file.txt -algorithm incremental -edge_selection 2 -initialization 1b
     Algorithm: incremental_edge_selection_2_initialization_1b
     Area: 100211210
     ratio: 0.263031
     Construction time: 5

•   ./to_polygon -i uniform-0000200-1.instance -o results_file.txt -algorithm incremental -edge_selection 2 -initialization 2a
     Algorithm: incremental_edge_selection_2_initialization_2a
     Area: 37363878
     ratio: 0.278152
     Construction time: 8

•   ./to_polygon -i protein-0020000.instance -o results_file.txt -algorithm incremental -edge_selection 3 -initialization 2b
     Algorithm: incremental_edge_selection_3_initialization_2b
     Area: 75886700638
     ratio: 0.685471
     Construction time: 151085

•   ./to_polygon -i paris-0000200.instance -o results_file.txt -algorithm incremental -edge_selection 3 -initialization 2b
     Algorithm: incremental_edge_selection_3_initialization_1a
     Area: 197092256
     ratio: 0.724213
     Construction time: 5

•   ./to_polygon -i euro-night-0100000.instance -o results_file.txt -algorithm incremental -edge_selection 1 -initialization 1b
     Algorithm: incremental_edge_selection_1_initialization_1b
     Area: 2332581614
     ratio: 0.407177
     Construction time: 4199035



    -- Για τον αλγόριθμο με βάση το ΚΠ --

 •   ./to_polygon -i stars-0000200.instance -o results_file.txt -algorithm convex_hull -edge_selection 1
     Algorithm: convex_hull_edge_selection_1
     Area: 159720393928
     ratio: 0.48022
     Construction time: 2020

 •   ./to_polygon -i london-0000100.instance  -o results_file.txt -algorithm convex_hull -edge_selection 2
     Algorithm: convex_hull_edge_selection_2
     Area: 91436178
     ratio: 0.239999
     Construction time: 192

 •   ./to_polygon -i euro-night-0000050.instance euro -o results_file.txt -algorithm convex_hull -edge_selection 2
     Algorithm: convex_hull_edge_selection_2
     Area: 18837314
     ratio: 0.287966
     Construction time: 13


•    ./to_polygon -i paris-0000200.instance -o results_file.txt -algorithm convex_hull -edge_selection 3
     Algorithm: convex_hull_edge_selection_3
     Area: 242592482
     ratio: 0.891403
     Construction time: 1742


    -i το input file (.instance με κατάλληλη δομή με tabs ανάμεσα)
    -o το output file που θα εκτυπωθούν τα δεδομένα όπως αναφέρει η εκφώνηση
    - algorithm για την εκτέλεση συγκεκριμένου αλγορίθμου (incremental or convex_hull)
    - edge_selection για την επιλογή κατάλληλης ακμής/σημείου (random ή πρόσθεση ελάσχιστου/μέγιστου εμβαδού)
    - initialization για τον αυξητικό αλγόριθμο για την επιλογή του sorting των σημείων (1a , 1b, 2a, 2b)


    -- Παρατηρήσεις --
• Η κεντρική διαφορά μεταξύ των αλγορίθμων είναι ότι ο αυξητικός περιορίζει ανά επανάληψη τον έλεγχο των ακμών που θα προστεθούν και μάλιστα επιλέγει ένα σημείο με βάση το αρχικό sorting που έχει γίνει.Σε αντίθεση με τον αλγόριθμο με βάση το ΚΠ ο οποιος ελέγχει σε καθε επανάληψη όλες τις ακμές του πολυγώνου και όλα τα υπολοιπόμενα σημεία για την εύρεση ορατότητας και εγγύτητας.Τα παραπάνω τεκμαίρωνται και από τον χρόνο εκτέλεσης του προγράμματος χρησιμοποιώντας τους παραπάνω αλγορίθμους.

• Επιπλέον, παρατηρείται οτι κατά την εκτέλεση εύρεσης τοπικού μέγιστου και τοπικού ελάχιστου πολύγωνου ο αλγόριθμος με βάση το ΚΠ αν και πιο αργός βρίσκει πολύγωνα με μικρότερο/μεγαλύτερο τοπικό ελάχιστο/μέγιστο εμβαδό σε σχέση με τον αυξητικό αλγόριθμο.

                    -- Εργασια 2 --
-- localSearch.hpp - localSearch.cpp --
Περιέχει κλάση που χρησιοποιείται για το local search και εχει μεθόδους που χρησιμοποιούνται για την εκτέλεση του αλγορίθμου και την επιστροφή του ratio και του area σημαντικό είναι οτι για προστέθηκαν οι μέθοδοι applyKPathRemoval και applyBlueRemoval για την τροποποίηση του τμήματος του μονοπατιού και της ακμής που αφαιρείται αντίστοιχα στην κλαση Polygonization.Αξίζει να σημειωθεί επιπλέον ότι για να μετράμε σωστά το εμβαδόν της προστηθέμενης και της αφαιρούμενης επιφάνειας χρησιμοποιούμε το CGAL::ON_BOUNDED_SIDE.Υπάρχουν επιπλέον ιδιοτητικές μέθοδοι της κλάσης για εύρεση των αλλαγών και εφαρμογή τους.

-- simulatedAnnealing.hpp  simulatedAnnealing.cpp spatialSubdivision.cpp--
Περιέχουν την κλάση που χρησιμοποιείται για simulated anealing.Μέσα στην κλάση υπάρχουν μέθοδοι για την εκτελεση του αλογριθμου και την ευρεση του εμβαδουκαι του ratio,καθως επισης και ιδιοτηκες κλασεις για την εκτελεση global και local αλλαγων,για ελεγχο εαν οι αλλαγες ειναι εγυρες ,αντικατασταη των ακμων και ευρεση αλλαγων.Επιπλεον περιεχονται μεθοδοι για αρχικοποιηση του kdtree και ευρεση των intersected ακμων μεσα σε μια αναζητουμενη περιοχη.Τελος υπαρχουν μεθοδοι για τμηματοποιση των σημειων εφαρμογοι των αλογριθμων σε αυτα και συνχονευση των υποπολυγωνων.


