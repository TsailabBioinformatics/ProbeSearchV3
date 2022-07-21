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
    "717V5": "data/717V5/g717v5_h1h2.fa",
    "PtrichocarpaV3.1": "data/PtrichocarpaV3.1/Ptrichocarpa_444_v3.0.fa",
    "PtrichocarpaV4.0": "data/PtrichocarpaV4.0/Ptrichocarpa_533_v4.0.fa",
    "DeltoidesWV94": "data/DeltoidesWV94/PdeltoidesWV94_445_v2.0.fa",
    "sPta717V1" : "data/sPta717V1/Pta717s_v1.1.fa",
    "sPta717tV1" : "data/sPta717tV1/sPta717tremula.fasta",
    "sPta717aV1" : "data/sPta717aV1/sPta717alba.fasta"
};

/**
 * listens for put request: calls bowtie2 child process and sends back alignment visualization
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
                // parse standard output (sam file)
                console.log(stdout);
                stdout = parse(String(req.body.read), String(req.body.db), stdout);
                res.send(stdout);
        }); 
    } else {
        let ts = Date.now();
        execSync('echo ">seq\n' + String(req.body.read) + '" > ' + ts + "input.fa"); // make input
        exec('./batmis/scripts/batmap -q ./' + ts + 'input.fa -g indices/' + String(req.body.db) + '/' + String(req.body.db) + ' -n' + String(req.body.mismatches) + ' -m50 -o ' + ts + 'output.txt', 
            (error, stdout, stderr) => {
                if (error) {
                    console.error(`exec error: ${error}`);
                    console.error(`standard error: ${stderr}`);
                    execSync('rm *bin');
                }
                // remove header from sam file
                exec('grep -v ^@ ' + ts + 'output.txt', (error, stdout, stderror) => {
                    console.log(stdout);
                    // remove input and output files 
                    execSync('rm ' + ts + 'input.fa ' + ts + 'output.txt');
                    // parse standard output (sam file)
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
    // sam = sort_sam(sam); // sorts SAM chronologically 
    for (var i = 0; i < sam.length - 1; i++) {
        target = sam[i].split("\t");
        var strand = "+";
        var complement = false;
        if (target[2] != "*") { // valid target
            if (target[1] === "16" || target[1] === "272") {
                strand = "-";
                complement = true;
            }
            res += "\t" + target[2] + " : " + target[3] + " (" + strand + ")\n"; 
            /* parse SAM file - get CIGAR, read, and reference sequence */
            var gene = get_gene(target[2], parseInt(target[3]), (parseInt(target[3]) + parseInt(sequence.length - 1)), db);
            res += "\tgene: " + gene + "\n";
            cigar = target[5];
            read = target[9];
            reference = String(execSync('samtools faidx ' +  db_dictionary[db] + ' ' + target[2] + ":" + target[3] 
                                        + "-" + (parseInt(target[3]) + sequence.length - 1) + ' | sed 1d')).toUpperCase();
            reference = reference.replace("\n", "");
            if (strand === "-") {
                read = read.split("").reverse().join("").replace("\n", ""); /* TODO: figure why we need to use 'replace()' here */
                reference = reference.split("").reverse().join("").replace("\n", "");
            }
            if (complement === true) {
                read = comp(read);
                reference = comp(reference);
            }
            res += illustrate(cigar, read, reference); 
        } else {
            res += "no targets found\n\n";
            break;
        } // if
    } // for
    res = sort_illustration(res);
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
    var illustration = "\t\tQ   " + read.substring(0, min) + "\n \t\t    "; // string illustrating alignment
    for (var i = 0; i < min; i++) {
        if (read.substring(i, i + 1) == reference.substring(i, i + 1)) {
            illustration += "|";
        } else {
            illustration += " ";
            mismatches++;
        } // if
    } // for
    illustration += "\n\t\tT   " + reference.substring(0, min) + "\n\t";
    illustration += "total mismatches: " + mismatches + "\n\n";
    return illustration;

} // illustrate


/**
 * takes in a sequence and returns its base-pair complement 
 * @param str input read or reference 
 * @returns complement
 */
function comp(str) {
    var res = "";
    for (var i = 0; i < str.length; i++) {
        if (str.charAt(i) === "A") { res += "T" }
        else if (str.charAt(i) === "T") { res += "A"}
        else if (str.charAt(i) === "C") { res += "G"}
        else { res += "C"}
    }
    return res;
} // comp


/**
 * sorts SAM files chronologically by chromosome coordinate
 * selection sort
 * @param sam SAM file (array of rows)
 * @returns sorted SAM file by coordinates 
 */
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
} // sort_sam


/**
 * takes in the illustration created from standard output,
 * and sorts its components chronologically by mismatch, then chromosome.
 * @param illustration 
 * @returns sorted illustration
 */
function sort_illustration(illustration) {
    var targets = illustration.split("\n\n");
    var sorted = [];
    var min = 0;
    
    sorted[0] = new Array(0);
    for (var i = 0; i < targets.length - 1; i++) {
        var target_info = targets[i].split("\n");
        var chromosome = target_info[0].split(" ")[0].replace(/\D/g, ""); // grabs chromosome coordinate 
        var mismatch = target_info[5].substr(-2).trim();
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
 * gets the gene model name from the gff3 file
 * @param chrom chromosome 
 * @param coord starting coordinate 
 * @returns gene model name 
 */
function get_gene(chrom, start_coord, end_coord, db) {
    try { 
        var chrom_table = String(execSync('grep ' + chrom + ' ./data/' + db + '/' + db + '.gene.table')).split("\n");
        var min = Math.abs(start_coord - parseInt(chrom_table[0].split("\t")[1]));
        var hit = 0;

        for (var i = 0; i < chrom_table.length; i++) {
            var chrom_row = chrom_table[i].split("\t");
            /* find nearest coordinate */
            if (Math.abs(start_coord - parseInt(chrom_row[1])) < min) { 
                min = Math.abs(start_coord - parseInt(chrom_row[1]));
                hit = i;
            }
            if (Math.abs(start_coord - parseInt(chrom_row[2])) < min) { 
                min = Math.abs(start_coord - parseInt(chrom_row[2]));
                hit = i;
            }
            if (Math.abs(end_coord - parseInt(chrom_row[1])) < min) {
                min = Math.abs(start_coord - parseInt(chrom_row[1]));
                hit = i;
            }
            if (Math.abs(end_coord - parseInt(chrom_row[2])) < min) {
                min = Math.abs(start_coord - parseInt(chrom_row[1]));
                hit = i;
            }
        } // for 

        /* determine where target hits */
        var hit_row = chrom_table[hit].split("\t");
        var gene = "";
        if (end_coord < parseInt(hit_row[1])) { // intergenic 
            gene = "intergenic";
            /* check if start pokes into previous hit - partial hit
            if (start_coord <= parseInt(chrom_table[hit - 1].split("\t")[2])) {
                gene = chrom_table[hit - 1].split("\t")[3];
            } */
        } else if (start_coord > parseInt(hit_row[2])) { // intergenic 
            gene = "intergenic";
            /* check if end pokes into next hit - partial hit
            if (end_coord >= parseInt(chrom_table[hit + 1].split("\t")[1])) {
                gene = chrom_table[hit + 1].split("\t")[3];
            } */
        } else if (start_coord >= parseInt(hit_row[1]) && end_coord <= parseInt(hit_row[2])) { // perfect hit
            gene = chrom_table[hit].split("\t")[3];
        } else { 
            gene = chrom_table[hit].split("\t")[3];    
        }
        return gene;
    } catch (e) {
        var gene = "N/A";
        return gene;
    }

} // get_gene





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
