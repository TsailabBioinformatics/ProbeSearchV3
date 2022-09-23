<script>
import SequenceForm from './components/SequenceForm.vue'
import AlignmentResult from './components/AlignmentResult.vue'
import FocusedResult from './components/FocusedResult.vue'
import Deck from './components/Deck.vue'
import axios from 'axios'
export default {
  name: "Home",
  components: {
    SequenceForm,
    AlignmentResult,
    FocusedResult,
    Deck
  },
  data() {
    return {
      count: 1,      // number of alignment results 
      results: [],   // result data
      headers: [],   // result headers
      show: 0,       // result to show
      fullscreen: false,  // full-screen mode
      deck: false         // card deck mode
    }
  },
  methods: {
    async run() { 
      const payload = {
              read: this.$refs.form.$data.read,
                db: "sPta717V2.0",
        mismatches: this.$refs.form.$data.mismatches
      };
      const res = await axios.put('/', payload);
      /* make header */
      this.headers[0] = ("input read:\t\t" + this.$refs.form.$data.read + "\n" + 
                          "read length:\t" + this.$refs.form.$data.read.length + "\n" + 
                          "database:\t\tsPta717V2.0\n");
      /* populate alignment results */
      this.results[0] = res.data; 
      this.show = 1; // show the first result
      this.$refs.form.$data.loading = false; // stop loading
    },
    switch_screen(screen) {
      const page = document.getElementById('page');
      page.classList.add('focused');
      const main = document.getElementById('main');
      main.classList.add("blur");
      if (screen === "fullscreen") { this.fullscreen = true; }
      else { this.deck = true; }
      window.pageYOffset || document.documentElement.scrollTop || document.body.scrollTop;
    },
    minimize_screen() {
      const main = document.getElementById('main');
      main.classList.remove("blur");
      const page = document.getElementById('page');
      page.classList.remove('focused');
      this.fullscreen = false;
      this.deck = false;
      window.pageYOffset || document.documentElement.scrollTop || document.body.scrollTop;
    }
  }
}
</script>


<template>
  <div id="page" style="width: 100%; margin: 0; display: flex; flex-direction: column; align-items: center; justify-content: center; ">
    <div id="main" class="main">

      <SequenceForm ref="form" v-on:valid="run()" />
      <div v-for="n in this.count" :key=n v-show="n === this.show">
        <AlignmentResult :id=n :data="results[n - 1]" :total="this.count" :header="headers[n - 1]"
                         v-on:swap-left="this.show--" v-on:swap-right="this.show++" 
                         v-on:fullscreen="switch_screen(`fullscreen`)" v-on:deck="switch_screen(`deck`)" /> 
      </div>
    </div>
    <!-- background components, toggled by AlignmentResult -->
    <FocusedResult v-if="this.fullscreen == true" :header="headers[this.show - 1]" :data="results[this.show - 1]"
                   v-on:back="minimize_screen()" /> 
    <Deck v-if="this.deck == true" :count="this.count" :headers="headers" :results="results" 
          v-on:back="minimize_screen()" />
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
    -ms-overflow-style: none; 
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
    scrollbar-width: auto; 
}
.main {
    place-items: flex-start; 
    flex-wrap: wrap;
} }

</style>