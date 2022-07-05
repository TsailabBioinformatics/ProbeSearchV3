<!-- 
  ProbeSearchV3!
    TODO
      - Netlify Large Media -> git LFS || deploy via in-house server
      - download results component 
-->

<script>
import SequenceForm from './components/SequenceForm.vue'
import AlignmentResult from './components/AlignmentResult.vue'
import FocusedResult from './components/FocusedResult.vue'
import axios from 'axios'
export default {
  name: "Home",
  components: {
    SequenceForm,
    AlignmentResult,
    FocusedResult
  },
  data() {
    return {
      count: 1,      // number of alignment results 
      results: [`input read:	AGCCACATGGCTCAATGGGAGAGTGCTCGACTAGAGGCCGAAGCCA
read length:	46
database:	717V5

target 1 - T08 : 5390892 (+)
    Q:	AGCCACATGGCTCAATGGGAGAGTGCTCGACTAGAGGCCGAAGCCA
    	||||||||||||||||||||||||||||||||||||||||||||||
    T:	AGCCACATGGCTCAATGGGAGAGTGCTCGACTAGAGGCCGAAGCCA
total mismatches: 0

target 2 - A08 : 5554895 (+)
    Q:	AGCCACATGGCTCAATGGGAGAGTGCTCGACTAGAGGCCGAAGCCA
    	||||||||||||||||||||||||||||||||||||||||||||||
    T:	AGCCACATGGCTCAATGGGAGAGTGCTCGACTAGAGGCCGAAGCCA
total mismatches: 0
`],   // array of result data
      show: 1,        // result to show
      focused: false
    }
  },
  methods: {
    async run() { 
      this.count = this.$refs.form.$data.db.length;
      for (var i = 0; i < this.count; i++) {
        const payload = {
          read: this.$refs.form.$data.read,
            db: this.$refs.form.$data.db[`${i}`]
        };
        const res = await axios.put('/', payload);
        this.results[i] = res.data; // populate alignment results
      } // for
      this.show = 1; // show the first result
    },
    focus_view() {
      const page = document.getElementById('page');
      page.classList.add('focused');
      const main = document.getElementById('main');
      main.classList.add("blur");
      this.focused = true;
      window.scrollTo(0, 0);
    },
    unfocus_view() {
      const main = document.getElementById('main');
      main.classList.remove("blur");
      const page = document.getElementById('page');
      page.classList.remove('focused');
      this.focused = false;
      window.scrollTo(0, 0);
    }
  }
}
</script>


<template>
  <div id="page" style="width: 100%; margin: 0; display: flex; flex-direction: column; align-items: center; justify-content: center; ">
    <div id="main" class="main">

      <SequenceForm ref="form" v-on:valid="run()" />
      <div v-for="n in this.count" :key=n v-show="n === this.show">
        <AlignmentResult :id=n :data="results[n - 1]" :total="this.count"
                         v-on:swap-left="this.show--" v-on:swap-right="this.show++" v-on:expand="focus_view()" /> 
      </div>
    </div>
    <FocusedResult v-if="this.focused == true" :id="this.show" :data="results[this.show - 1]"
                   v-on:contract="unfocus_view()" /> 
  </div>
</template>

<style>
@import './assets/base.css';
#app {
  height: 100vh;
  width: 100vw;
  display: flex;
  align-items: start;
  overflow: auto;
  -ms-overflow-style: none; /* for Internet Explorer, Edge */
  scrollbar-width: none; 
}
.main {
  display: flex; 
  flex-direction: column;
  align-items: center; 
  justify-content: center;
}
.focused {
  display: grid;
  width: 100vw;
  min-height: 100vh;
}
.blur {
  display: none;
}
@media (min-width: 1024px) {
  body {
    display: flex;
    align-items: flex-start;
    justify-content: center;
  }
  span:hover {
    cursor: pointer;
  }
  #app {
    background-color: transparent;
    font-weight: normal;
    align-self: center;
    align-items: flex-start;
  }
 .main {
    place-items: flex-start; 
    flex-wrap: wrap;
  }
}
</style>





<!-- 
  CODE
      dynamically adding forms 

      for (var i = 0; i <= this.extra; i++) {

        if (i === 0) {
            const payload = {
              read: this.$refs[`${i}`].$data.read,
                db: this.$refs[`${i}`].$data.db
              };
            const res = await axios.put('/', payload);
            this.showres = true
            this.res += `Search Result ${i+1}\n-----------------\n`
            this.res += res.data + '\n';
        } else {
            console.log(this.$refs[`${i}`][0].$data.read);
            const payload = {
              read: this.$refs[`${i}`][0].$data.read,
              db: this.$refs[`${i}`][0].$data.db
              };
            const res = await axios.put('/', payload);
            this.res += `Search Result ${i+1}\n-----------------\n`

            this.res += res.data + '\n';
        } // if/else

      } // for

      <div v-for="i in extra" :key="i">
        <SequenceForm :ref="`${i}`" :mainform="false" v-on:subtract-form="this.extra--"/>
      </div>
      v-on:add-form="this.extra++"

      .add-form {
        color: darkgrey;
        position: absolute;
        right: 100px; 
        bottom: 100px;
        cursor: pointer;
      }
      .hidden-instructions {
        display: none;
      }
      .add-form:hover + .hidden-instructions {
        display: flex;
        padding: 5px;
        width: 110px;
        font-size: 10px;
        text-align: center;
        color: rgba(60, 60, 60, 0.66);
        background-color: rgba(247, 247, 247, 0.90);
        border: 1px solid hsla(160, 100%, 37%, 1);
        position: absolute;
        right: -6.5px;
        bottom: -25px;
      }


    .add-form {
      color: darkgrey;
      position: absolute;
      right: 10px; 
      bottom: 0;
    }
    .add-form:hover {
      cursor: pointer;
    }
    .hidden-instructions {
      color: red;
    }
    .add-form:hover + .hidden-instructions {
      display: flex;
      color: red;
    }
    */
    /*
    @media (hover: hover) {
      a:hover {
        background-color: hsla(160, 100%, 37%, 0.2);
      }
    }
-->
