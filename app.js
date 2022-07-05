const express = require('express');
const path = require('path');
const { exec } = require('child_process');
const { execSync } = require('child_process');
const { stdout } = require('process');

const public = __dirname + '/vue/dist/';
const app = express();

app.use(express.json());
app.use(express.static(public))

app.get('/', function (req, res) {
    res.sendFile(path.join(public + "index.html"));
});

app.listen({ port: 8080 }, async () => {
    console.log('Server up @ http://localhost:8080/ !')
})

/* equates user-selected db name to its path in the directory */
var db_dictionary = {
    "717V5": "data/g717v5_h1h2.fa",
    "PtrichocarpaV3.1": "data/Ptrichocarpa_444_v3.0.fa",
    "PtrichocarpaV4.0": "data/Ptrichocarpa_533_v4.0.fa",
    "DeltoidesWV94": "data/PdeltoidesWV94_445_v2.0.fa",
    "sPta717V1" : "data/Pta717s_v1.1.fa",
    "sPta717tV1" : "data/sPta717tremula.fasta",
    "sPta717aV1" : "data/sPta717alba.fasta"
};

/**
 * listens for put request: calls bowtie2 child process and sends back alignment visualization
 */
app.put('/',  async (req, res) => {
    /**
     * calls bowtie2 
     */
    exec('./bowtie/bowtie2 -x indices/' + String(req.body.db) + ' -k 30 -c ' + String(req.body.read) + '  --end-to-end --no-hd' , (error, stdout, stderr) => {
        if (error) {
          console.error(`exec error: ${error}`);
        }
        // parse standard output for nice visualization
        console.log(stdout);
        stdout = parse(String(req.body.read), String(req.body.db), stdout);
        res.send(stdout);
    });
})

/**
 * parses the SAM file input and gathers reference sequence.
 * @param sequence input read 
 * @param db database
 * @param sam SAM file  
 */
function parse(sequence, db, sam) {
    res = "input read:\t" + sequence + "\nread length:\t" + sequence.length + "\ndatabase:\t" + String(db) + "\n\n"; // instatiate string
    sam = sam.split("\n");
    for (var i = 0; i < sam.length - 1; i++) {
        target = sam[i].split("\t");
        var strand = "+";
        if (target[2] != "*") { // valid target
            if (target[1] === "16") {
                strand = "-";
            }
            res += "target " + String(i + 1) + " - " + target[2] + " : " + target[3] + " (" + strand + ")\n"; 
            /* parse SAM file - get CIGAR, read, and reference sequence */
            cigar = target[5];
            read = target[9];
            reference = String(execSync('samtools faidx ' +  db_dictionary[db] + ' ' + target[2] + ":" + target[3] 
                                        + "-" + (parseInt(target[3]) + sequence.length - 1) + ' | sed 1d'));
            reference = reference.replace("\n", "");
            if (strand === "-") {
                read = read.split("").reverse().join("").replace("\n", ""); /* TODO: figure why we need to use 'replace()' here */
                reference = reference.split("").reverse().join("").replace("\n", "");
            }
            res += illustrate(cigar, read, reference); 
        } else {
            res += "no targets found\n\n";
            break;
        } // if
    } // for
    console.log(res)
    return res;
} // parse 

/**
 * illustrates the alignment between the read & reference.
 * @param cigar CIGAR string of SAM
 * @param read sequence of read 
 * @param reference sequence of reference 
 * @returns illustration of alignment
 */
function illustrate(cigar, read, reference) {
    
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
            var space = " ".repeat(nucleotides);
            reference = reference.substring(0, ptr) + space + reference.substring(ptr, reference.length);
            ptr += nucleotides;
        } else if (letter === "D") {
            var space = " ".repeat(nucleotides);
            read = read.substring(0, ptr) + space + read.substring(ptr, read.length);
            ptr += nucleotides;
        }
        if (j != cigar.length - 1) {
            i = j;
        } // if
    } // for
    
    /* build illustration  */
    var min = Math.min(read.length, reference.length);
    var illustration = "    Q:\t" + read.substring(0, min) + "\n    \t"; // string illustrating alignment
    for (var i = 0; i < min; i++) {
        if (read.substring(i, i + 1) == reference.substring(i, i + 1)) {
            illustration += "|";
        } else {
            illustration += " ";
            mismatches++;
        } // if
    } // for
    illustration += "\n    T:\t" + reference.substring(0, min) + "\n";
    illustration += "total mismatches: " + mismatches + "\n\n";
    return illustration;

} // illustrate

function reverse(str) {
    res = str.split("").reverse().join("")
}













/*
notes
`````
insertion = extra nucleotide in read - add space in reference 
            eg, 3M2I3M:
            read: A T G C A T G C
             ref: A T G     T G C
                

deletion = extra nucleotide in reference - add space in read 
            eg, 3M2D3M:
            read: A T G     T G C
             ref: A T G G C T G C


// other possible CIGAR letters //
cigar.substring(i, cigar.length).indexOf("D"), cigar.substring(i, cigar.length).indexOf("N"),
cigar.substring(i, cigar.length).indexOf("S"), cigar.substring(i, cigar.length).indexOf("H"),
cigar.substring(i, cigar.length).indexOf("P")); 

*/


/*
code
````
*/
