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
Περιέχει κλάση που χρησιοποιείται για το local search και εχει μεθόδους που χρησιμοποιούνται για την εκτέλεση του αλγορίθμου(start) και την επιστροφή του ratio και του area(getOptimisedArea,getOptimisedRatio) σημαντικό είναι οτι για προστέθηκαν οι μέθοδοι applyKPathRemoval και applyBlueRemoval για την τροποποίηση του τμήματος του μονοπατιού και της ακμής που αφαιρείται αντίστοιχα στην κλαση Polygonization.Αξίζει να σημειωθεί επιπλέον ότι για να μετράμε σωστά το εμβαδόν της προστηθέμενης και της αφαιρούμενης επιφάνειας χρησιμοποιούμε το CGAL::ON_BOUNDED_SIDE.Υπάρχουν επιπλέον ιδιοτηκές μέθοδοι της κλάσης για εύρεση των αλλαγών και εφαρμογή τους(findChanges,applyChanges).

Εάν επιλέγει απο τον χρήστη να γινει local search τοτε εκτελείται ο local search σε ολο το πολύγωνο και αναζητούνται ολα τα μονοπάτια σημείων μέσα στο πολύγωνο που βρίσκονται στις ακμές του και αντικαθήστανται με άλλη ακμή του πολυγώνου διατειρόντας παράλληλα την απλότητα του πολυγώνου.

-- simulatedAnnealing.hpp  simulatedAnnealing.cpp spatialSubdivision.cpp--
Περιέχουν την κλάση που χρησιμοποιείται για simulated anealing.Μέσα στην κλάση υπάρχουν μέθοδοι για την εκτέλεση του αλογoρίθμου(startAnnealing ή startSubdivision) και την εύρεση του εμβαδού και του ratio(getOptimisedArea,getOptimisedRatio),καθώς επίσης και ιδιοτηκές κλάσεις για την εκτέλεση global και local αλλαγών(localTransition,globalTransition),για έλεγχο εάν οι αλλαγές είναι έγκυρες(validityCheck),αντικατάσταση των ακμών και εύρεση αλλαγών(replace,findGlobalChanges).Επιπλέον περιέχονται μέθοδοι για αρχικοποίηση του kdtree και εύρεση των intersected ακμών μέσα σε μια αναζητούμενη περιοχή(KdTreeInit,validityCheck).Τέλος υπάρχουν μέθοδοι για τμηματοποίηση των σημείων εφαρμογή των αλογορίθμων σε αυτά και συνχονευση των υποπολυγώνων(createSubsetPoints,mergePolygons).

Εάν επιλεγεί απο τον χρήστη το local anealing εφαρμόζονται local transitions σε όλο το πολύγωνο οπου επιλέγεται ένα τυχαίο σημείο και αλλάζει θέση με το επομένο του,η εγκυρότητα του πολυγώνου ελέγχεται με kd-tree.

Εάν επιλεγεί απο τον χρήστη το global anealing εφαρμόζονται local transitions σε όλο το πολύγωνο οπου επιλέγονται δύο σημεία ενόνωνται μεταξύ τους και στην συνέχεια διαγράφεται η ακμή του σημείου που ενώνεται με το πρώτο και στην συνέχεια προσθέτουμε νεο segment σε αυτη τη περιοχή,η εγκυρότητα του πολυγώνου και εδώ ελέγχεται με kd-tree.

Σε περίπτωση που επιλέγει απο τον χρήστη το subdision διαμερίζουμε τα σημεία σε k υποσύνολα με βάση τον τύπο k = ceil(n-1/m-1) το m ορίζεται από τον χρήστη.Το τελευταίο segment καθε πολυγώνου εκτός απο την τελευταία ομάδα απαιτείται να είναι γνησίως μονότονο,ενώ το πρώτο εκτός της πρώτης ομάδας απαιτείται αν είναι γνησιώς φθήνων.Εφαρμογή στην συνέχεια ενός εκ των αρχικών αλγορίθμων πολυγωνοποίησης(convex hull ή incremental) σε κάθε ομάδα σημείων,αποφυγή επεξεργασίας των άκρων των ομάδων καθώς μπορεί να προκληθεί σφάλμα κατά την συνχώνευση των υποομάδων.Εάν κατά την διαδικασία εκτέλεσης του incrental o αλγόριθμος αποτύχει να κρατήσει το marked segment τότε επιλέγουμε να εκτελεστεί ο convex hull για την συγκεκριμένη ομάδα.Στην συνέχεια εκτελούμε global transitions για κάθε υποπολύγωνο διατειρόντας τα ευθύγραμμα τμήματα συγχώνευσης.Για την συγχώνευσευση όλων των πολυγώνων ξεκινάμε από το lower hull ολων των πολυγώνων και ενώνουμε την κάθε ομάδα και μόλις φτάσουμε στο τέλος ενόνουμε τα uper hull των ομάδων πηγαίνοντας στην την φορά απο το τελός στην άρχη,στο τελευταίο στάδιο του merging γίνεται η διαγραφή των κοινών ακμών κάθε ομάδων με κοινό σημείο και η ένωση του κάτω μέρους του τριγώνου που προκύπτει.Τλος εφαρμόζουμε local transitions στο πολύγωνο που προκύπτει.


Παρατηρήσεις:
Χάρη στο subdivision πετυχαίνουμε καλύτερους χρόνους κατασκευής σε σχέση με την πρώτη εργασία,έχουμε επίσης κάνει optimise στην isSimple που ελέγχει εάν το πολύγωνο είναι απλό.
Κατά την εκτέλεση του simulated anealing μπορεί να προκύψει πολύγωνο με μικρότερο εμβαδόν εάν δοθεί maximazation επιλογή απο τον χρήστη ή και το αντίστροφο αυτό μπορεί να αποφευχθεί είτε αυξάνοντας το L(που και πάλι δεν είναι εγγυημένο οτι θα έρθουν καλύτερα αποτελέσματα) ή κρατόντας το στιγμειότυπο του πολυγώνου κάθε φορά που έχει αυτή την στιγμή το μεγαλύτερο δυνατό εμβαδόν.
Ο local search είναι πιό γρήγορος για λίγα σημεία ομως για πολλά (πάνω απο 200) αργεί αρκετά να δώσει απαντήσεις ειδικά για μικρο threshold, και μεγάλο K.

Ο simulated annealing κάνει καλύτερες προσεγγίσεις στο εμβαδόν από την άλλη και είναι ικανοποιητικός ως προς τον χρόνο ειδικά εαν επιλεγεί το subdivision.Τα local transitions έχουν γρήγορα αποτελέσματα αλλά όχι τόσο optimized,τα global transitions είναι αρκετά πιο αργά σε σχέση με τον local αλλα τα αποτελέματα του είναι αρκετά πιο optimal.Τέλος ο spatial subdivision συνδιάζοντας global transition σε μικρά πολύγωνα και local transition στο τελικό πολύγωνο παρατειρούμε ένα αρκετά ικανοποιητικό αποτέλεσμα ως προς τον χρόνο και ως προς την προσέγγιση. 

Μεταγλώττιση:

• cgal_create_CMakeLists -s optimal_polygon
• cmake -DCMAKE_BUILD_TYPE=Release .
• make

Εκτέλεση
./optimal_polygon -i αρχειο_εισοδου -o αρχειο_εξοδου -algorithm αλγοριθμος_optimization(local_search or simulated_annealing) -L αριθμός -max ή -min -annealing μεθοδος_annealing(local ή global ή subdivision) -algorithm_initial αρχικός_αλγόριθμος(incremental ή convex_hull) -initialization μονο_για_τον_incremental (1a ή 1b ή 2a ή 2b. Ο subdivision μόνο με 1a και 1b) -m αριθμός -edge_selection (1 random , 2 min, 3 max, αν δεν δωθεί τότε εκτελεί οτι κάνει το optimisation)


Παραδείγματα εκτέλεσης

./optimal_polygon -i uniform-0000030-1.instance -o output.txt -algorithm local_search -L 3 -min -threshold 2.0 -algorithm_initial convex_hull
Algorithm: local_search_min
Area_initial: 802586
Area: 550732
ratio_initial: 0.377536
ratio: 0.259064
Construction time: 44


./optimal_polygon -i uniform-0000500-1.instance -o output.txt -algorithm simulated_annealing -L 100000 -max  -annealing local -algorithm_initial incremental -initialization 1a
Algorithm: simulated_annealing_max
Area_initial: 606856524
Area: 644944060
ratio_initial: 0.698769
ratio: 0.742625
Construction time: 3020


./optimal_polygon -i euro-night-0000300.instance -o output.txt -algorithm simulated_annealing -L 6000 -min -annealing global -algorithm_initial incremental -initialization 1a
Algorithm: simulated_annealing_min
Area_initial: 18620510
Area: 12292822
ratio_initial: 0.229772
ratio: 0.15169
Construction time: 18316


./optimal_polygon -i uniform-0000100-1.instance -o output.txt -algorithm local_search -L 5 -max -threshold 2.0 -algorithm_initial incremental -initialization 1a
Algorithm: local_search_max
Area_initial: 23777634
Area: 26522038
ratio_initial: 0.713752
ratio: 0.796133
Construction time: 4474


./optimal_polygon -i uniform-0000100-1.instance -o output.txt -algorithm simulated_annealing -L 10000 -max -annealing global -algorithm_initial incremental -initialization 1a
Algorithm: simulated_annealing_max
Area_initial: 23777634
Area: 26191476
ratio_initial: 0.713752
ratio: 0.78621
Construction time: 4448


./optimal_polygon -i uniform-0000100-1.instance -o output.txt -algorithm simulated_annealing -L 10000 -max -annealing global -algorithm_initial incremental -initialization 1a -edge_selection 1
Algorithm: simulated_annealing_max
Area_initial: 18128506
Area: 25923254
ratio_initial: 0.544178
ratio: 0.778159
Construction time: 4077


./optimal_polygon -i uniform-0000500-1.instance -o output.txt -algorithm simulated_annealing -L 6000 -max -annealing subdivision -algorithm_initial convex_hull -m 20
Algorithm: simulated_annealing_max
Area_initial: 559346838
Area: 605628460
ratio_initial: 0.644063
ratio: 0.697355
Construction time: 8061


./optimal_polygon -i uniform-0000500-1.instance -o output.txt -algorithm simulated_annealing -L 6000 -max -annealing subdivision -algorithm_initial incremental -initialization 1a -m 20
Algorithm: simulated_annealing_max
Area_initial: 535270110
Area: 608276158
ratio_initial: 0.61634
ratio: 0.700403
Construction time: 7012

