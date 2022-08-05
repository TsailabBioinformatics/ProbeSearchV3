const express = require('express');
const path = require('path');
const { exec } = require('child_process');
const { execSync } = require('child_process');
const sqlite = require('better-sqlite3');

const public = __dirname + '/vue/717v2/dist/';
const app = express();

app.use(express.json());
app.use(express.static(public))

app.get('/', function (req, res) {
    res.sendFile(path.join(public + "index.html"));
});

app.listen({ port: 8888 }, async () => {
    console.log('server up @ http://localhost:3000/ !')
});

let db = new sqlite('data/sPta717V2.0/coverage.db');


/*********************************************************************************************************************/


/**
 * listens for put request: calls batmis or bowtie2 child process and sends back alignment visualization
 */
app.put('/',  async (req, res) => {
    /**
     * call bowtie2 or batmap
     */
    if (req.body.read.length > 50) {
        exec('./bowtie/bowtie2 -x indices/sPta717V2.0/sPta717V2.0 -k 30 -c ' + String(req.body.read) + '  --very-sensitive --no-hd' , 
            (error, stdout, stderr) => {
                if (error) {
                    console.error(`exec error: ${error}`);
                }
                // parse sam file
                console.log(stdout);
                stdout = parse(String(req.body.read), "sPta717V2.0", stdout);
                res.send(stdout);
        }); 
    } else {
        let ts = Date.now();
        // make fasta input
        execSync('echo ">seq\n' + String(req.body.read) + '" > ' + ts + "input.fa"); 
        exec('./batmis/scripts/batmap -q ./' + ts + 'input.fa -g indices/sPta717V2.0/sPta717V2.0 -n' + String(req.body.mismatches) + ' -m50 -o ' + ts + 'output.txt', 
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
                    stdout = parse(String(req.body.read), "sPta717V2.0", stdout);
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
            /* parse SAM file - get CIGAR, read, and reference sequence */
            var gene = get_gene(target[2], parseInt(target[3]), (parseInt(target[3]) + parseInt(sequence.length - 1)), db);
            res += "\tgene: " + gene + "\n";
            cigar = target[5];
            read = target[9];
            reference = String(execSync('samtools faidx data/sPta717V2.0/sPta717_v2.0.fa ' + target[2] + ":" + target[3] 
                                        + "-" + (parseInt(target[3]) + sequence.length - 1) + ' | sed 1d')).toUpperCase();
            reference = reference.replace("\n", "");
            if (strand === "-") {
                read = read.split("").reverse().join("").replace("\n", ""); 
                reference = reference.split("").reverse().join("").replace("\n", "");
            }
            if (complement === true) {
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
    var ptr = 0, mismatches = 0; 
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

        var letter = cigar.charAt(j), nucleotides = parseInt(cigar.substring(i, j));
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
        if (j != cigar.length - 1) { i = j; } 
    } // for
    
    /* build illustration  */
    var illustration = "\t\tQ   " + read.substring(0, read.length) + "\n \t\t    "; 
    for (var i = 0; i < read.length; i++) {
        if (read.substring(i, i + 1) == reference.substring(i, i + 1)) { illustration += "|"; } 
        else { illustration += " "; mismatches++; } 
    } // for
    return illustration += "\n\t\tT   " + reference.substring(0, read.length) 
                        + get_coverage(chrom, parseInt(pos), parseInt(read.length)) 
                        + "\n\t" + "total mismatches: " + mismatches + "\n\n";

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
        var chromosome = target_info[0].split(" ")[0].replace(/\D/g, ""); // grabs chromosome coordinate 
        var mismatch = target_info[(x + 5)].substr(-2).trim();
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
    } // for
    return res;
} // sort_illustration


/**
 * gets the gene model name from the gff3 file via `bedtools intersect`
 * @param chrom chromosome 
 * @param start_coord starting coordinate
 * @param end_coord end coordinate  
 * @returns gene model name 
 */
 function get_gene(chrom, start_coord, end_coord, db) {
    let ts = Date.now();
    execSync('echo "' + String(chrom) + '\t' + String(start_coord) + '\t' + String(end_coord) + '" > ' + ts + "query.bed") // make input file
    var intersect = String(execSync('bedtools intersect -a ' + ts + 'query.bed -b data/' + db + '/' + db + '.gene.table -wb -nonamecheck | head -n 1'));
    execSync('rm ' + ts + 'query.bed'); // delete input file 

    if (intersect === "") { return "intergenic" }
    else { return intersect.split("\t")[6].replace("\n", "") }
} // get_gene


/**
 * gets the coverage score from 
 * @param chrom chromosome 
 * @param pos position
 * @param length length of read
 * @returns coverage score
 */
function get_coverage(chrom, pos, length) {
    var coverage = "\n\t\tC   ";
    let sql = `SELECT * FROM ${chrom} WHERE pos BETWEEN ${pos} AND ${(pos + (length - 1))} LIMIT ${length}`;
    if (chrom.charAt(0) === "s") { sql = `SELECT * FROM Scaffold WHERE chrom='${chrom}' AND pos BETWEEN ${pos} AND ${(pos + (length - 1))} LIMIT ${length}`; } 
    var rows = db.prepare(sql).all();

    if (rows.length != 0) {
        for (const row of rows) { 
            var score = row.score;
            if (score == 0) { coverage += "0" }
            else if (score <= 10) { coverage += "A" }
            else if (score <= 50) { coverage += "B" }
            else if (score <= 100) { coverage += "C" }
            else if (score <= 500) { coverage += "D" }
            else { coverage += "E" }
        }
    } else { coverage += "0".repeat(length); }
    return coverage;
}