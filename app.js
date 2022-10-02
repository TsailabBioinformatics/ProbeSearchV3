const express = require('express');
const path = require('path');
const { exec } = require('child_process');
const { execSync } = require('child_process');
const public = __dirname + '/vue/dist/';
const app = express();

app.use(express.json());
app.use(express.static(public))

app.get('/', function (req, res) {
    res.sendFile(path.join(public + "index.html"));
});

app.listen({ port: 8090 }, async () => {
    console.log('server up @ http://localhost:8090/ !')
});


/* equates user-selected db name to its path in the directory */
var db_dictionary = {
    "717V5": "data/717V5/g717v5_h1h2.fa",
    "PtrichocarpaV3.1": "data/PtrichocarpaV3.1/Ptrichocarpa_444_v3.0.fa",
    "PtrichocarpaV4.0": "data/PtrichocarpaV4.0/Ptrichocarpa_533_v4.0.fa",
    "DeltoidesWV94": "data/DeltoidesWV94/PdeltoidesWV94_445_v2.0.fa",
    "sPta717V1" : "data/sPta717V1/Pta717s_v1.1.fa",
    "sPta717tV1" : "data/sPta717tV1/sPta717tremula.fasta",
    "sPta717aV1" : "data/sPta717aV1/sPta717alba.fasta",
    "sPta717V2.0" : "data/sPta717V2.0/sPta717_v2.0.fa"
};


/*********************************************************************************************************************/


/**
 * listens for put request: calls batmis or bowtie2 child process and sends back alignment visualization
 */
app.put('/',  async (req, res) => {
    /**
     * call bowtie2 or batmap
     */
    if (req.body.read.length > 50) {
        exec('./bowtie/bowtie2 -x indices/' + String(req.body.db) + '/' + String(req.body.db) + ' -k 30 -c ' + String(req.body.read) + '  --very-sensitive --no-hd' , 
            (error, stdout, stderr) => {
                if (error) {
                    console.error(`exec error: ${error}`);
                }
                // parse sam file
                console.log(stdout);
                stdout = parse(String(req.body.read), String(req.body.db), stdout);
                res.send(stdout);
        }); 
    } else {
        let ts = Date.now();
        // make fasta input
        execSync('echo ">seq\n' + String(req.body.read) + '" > ' + ts + "input.fa"); 
        exec('./batmis/scripts/batmap -q ./' + ts + 'input.fa -g indices/' + String(req.body.db) + '/' + String(req.body.db) + ' -n' + String(req.body.mismatches) + ' -m50 -o ' + ts + 'output.txt', 
            (error, stdout, stderr) => {
                if (error) {
                    console.error(`exec error: ${error}`);
                    execSync('rm *bin');
                }
                // remove header from sam file
                exec('grep -v ^@ ' + ts + 'output.txt', (error, stdout, stderror) => {
                    console.log(stdout);
                    // remove input and output files 
                    execSync('rm ' + ts + 'input.fa ' + ts + 'output.txt');
                    // parse sam file
                    stdout = parse(String(req.body.read), String(req.body.db), stdout);
                    res.send(stdout);
                });                   
        }); 
    }
})

/**
 * parses the SAM file input and gathers reference sequence.
 * @param sequence input read 
 * @param db database
 * @param sam SAM file  
 */
function parse(sequence, db, sam) {
    res = "" 
    sam = sam.split("\n");
    for (var i = 0; i < sam.length - 1; i++) {
        target = sam[i].split("\t");
        var strand = "+";
        var complement = false;
        if (target[2] != "*") { // valid target
            if (target[1] === "16" || target[1] === "272") { // reverse complement
                strand = "-";
                complement = true;
            }
            res += "\t" + target[2] + " : " + target[3] + " (" + strand + ")\n"; 
            var gene = get_gene(target[2], parseInt(target[3]), (parseInt(target[3]) + parseInt(sequence.length - 1)), db);
            res += "\tgene: " + gene + "\n";
            cigar = target[5];
            read = target[9];
            reference = String(execSync('samtools faidx ' +  db_dictionary[db] + ' ' + target[2] + ":" + target[3] 
                                        + "-" + (parseInt(target[3]) + sequence.length - 1) + ' | sed 1d')).toUpperCase();
            reference = reference.replace("\n", "");
            if (strand === "-") {
                read = read.split("").reverse().join("").replace("\n", ""); 
                reference = reference.split("").reverse().join("").replace("\n", "");
                read = comp(read);
                reference = comp(reference);
            }
            res += illustrate(cigar, read, reference, db, target[2], target[3]); 
        } else {
            res += "no targets found\n\n";
            break;
        } // if
    } // for
    res = sort_illustration(res, 1); 
    return res;
} // parse 

/**
 * illustrates the alignment between the read & reference.
 * @param cigar CIGAR string of SAM
 * @param read sequence of read 
 * @param reference sequence of reference 
 * @returns illustration of alignment
 */
function illustrate(cigar, read, reference, db, chrom, pos) {
    
    var ptr = 0; // alignment iterator
    var mismatches = 0; // number of mismatches 
    /* parse cigar  */
    for (var i = 0; i < cigar.length; i++) { 

        var j = cigar.length - 1;
        // check letters and grab index of closest letter
        if (cigar.substring(i, cigar.length).indexOf('M') != -1) { 
            j = cigar.substring(i, cigar.length).indexOf('M') + i;
        }
        if (cigar.substring(i, cigar.length).indexOf('I') != -1 && (cigar.substring(i, cigar.length).indexOf('I') + i) < j) { 
            j = cigar.substring(i, cigar.length).indexOf('I') + i;
        }
        if (cigar.substring(i, cigar.length).indexOf('D') != -1 && (cigar.substring(i, cigar.length).indexOf('D') + i) < j) { 
            j = cigar.substring(i, cigar.length).indexOf('D') + i;
        }
        var letter = cigar.charAt(j);
        var nucleotides = parseInt(cigar.substring(i, j));

        if (letter === 'M') {
            ptr += nucleotides;
        } else if (letter === "I") {
            var space = "-".repeat(nucleotides);
            reference = reference.substring(0, ptr) + space + reference.substring(ptr, reference.length);
            ptr += nucleotides;
        } else if (letter === "D") {
            var space = "-".repeat(nucleotides);
            read = read.substring(0, ptr) + space + read.substring(ptr, read.length);
            ptr += nucleotides;
        }
        if (j != cigar.length - 1) {
            i = j;
        } // if
    } // for
    
    /* build illustration  */
    var min = Math.min(read.length, reference.length);
    var illustration = "\t\tQ   " + read.substring(0, min) + "\n \t\t    "; 
    for (var i = 0; i < min; i++) {
        if (read.substring(i, i + 1) == reference.substring(i, i + 1)) {
            illustration += "|";
        } else {
            illustration += " ";
            mismatches++;
        } // if
    } // for
    return illustration += "\n\t\tT   " + reference.substring(0, min) + "\n\t" + "total mismatches: " + mismatches + "\n\n";

} // illustrate


/**
 * takes in a sequence and returns its base-pair complement 
 * @param str input read or reference 
 * @returns complement
 */
function comp(str) {
    var res = ""; 
    var complement = { "A": "T", "T": "A", "G": "C", "C": "G" }
    for (var i = 0; i < str.length; i++) {
        res += complement[str.charAt(i)]
    } // for
    return res;
} // comp


/**
 * takes in the illustration created from standard output,
 * and sorts its components chronologically by mismatch, then chromosome.
 * @param illustration 
 * @returns sorted illustration
 */
function sort_illustration(illustration, x) {
    var targets = illustration.split("\n\n");
    var sorted = [];
    var min = 0;

    sorted[0] = new Array(0);
    for (var i = 0; i < targets.length - 1; i++) {
        var target_info = targets[i].split("\n");
        var chromosome = target_info[0].split(" ")[0].replace(/\D/g, ""); 
        var mismatch = target_info[(x + 4)].substr(-2).trim();
        if (mismatch != min) {
            sorted[mismatch] = new Array(0);
            min = mismatch;
        }
        sorted[mismatch].push([chromosome, targets[i]]);
    } // for
    
    var res = "";
    var count = 0;
    for (var j = 0; j < sorted.length; j++) {
        if (sorted[j]) { 
            sorted[j] = sorted[j].sort(); 
            for (var k = 0; k < sorted[j].length; k++) {
                res += "hit " + (++count) + ")\n" + sorted[j][k][1] + "\n\n"
            }
        }
    }
    return res;
}


/**
 * gets the gene model name from the gff3 file via `bedtools intersect`
 * @param chrom chromosome 
 * @param coord starting coordinate 
 * @returns gene model name 
 */
function get_gene(chrom, start_coord, end_coord, db) {
    let ts = Date.now();
    execSync('echo "' + String(chrom) + '\t' + String(start_coord) + '\t' + String(end_coord) + '" > ' + ts + "query.bed") 

    var intersect = String(execSync('bedtools intersect -a ' + ts + 'query.bed -b data/' + db + '/' + db + '.gene.table -wb -nonamecheck | head -n 1'));
    execSync('rm ' + ts + 'query.bed');

    if (intersect === "") { return "intergenic" }
    else { 
        console.log("hit")
        return intersect.split("\t")[6].replace("\n", "") 
    }
} // get_gene



/**************************************************************************************************/



/*
 * sorts SAM files chronologically by chromosome coordinate
 * selection sort
 * @param sam SAM file (array of rows)
 * @returns sorted SAM file by coordinates 

 function sort_sam(sam) {
    var min = 0;
    for (var i = 0; i < sam.length - 1; i++) {
        min = i;
        for (var j = i + 1; j < sam.length - 1; j++) {
            var x = sam[j].split("\t")[2]; // represents the chromosome in a row of the SAM file
            var y = sam[min].split("\t")[2];
            if (parseInt(x.substring(x.length - 2, x.length)) < parseInt(y.substring(y.length - 2, y.length))) {
                min = j;
            }
        }
        var temp = sam[min];
        sam[min] = sam[i];
        sam[i] = temp; 
    }
    return sam;
} // sort_sam */