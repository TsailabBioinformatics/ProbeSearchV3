<script>
  import Card from './Card.vue'
  export default {
    name: "Deck",
    components: {
        Card
    },
    props: ['count', 'headers', 'results'],
    methods: {
      goback() {
        this.$emit("back")
      },
    },
    emits: ["back"]
  }
</script>


<template>
  <div style="width: 100%"> 
    <h3 class="goback" @click="goback()"> ← </h3>
    <p class="instruction"> go back </p>  
  </div>

  <div class="board">
        
        <div class="spread">
            <Card v-for="n in count" :header="headers[n - 1]" :result="results[n - 1]" />
        </div>
    
  </div>

</template>


<style scoped>

.board {
    width: 90vw;
    height: 100%;
    min-height: 90vh;
    background-color: transparent;
}
.spread {
    display: flex;
    align-items: flex-start;
    justify-content: flex-start;
}
.goback {
    position: absolute; 
    left: 10px; 
    top: 0;
}
.goback:hover {
    cursor: pointer;
}
.instruction {
    display: none;
}

@media (min-width: 1024px) {
.board {
    margin: 5vh 5vw 5vh 5vw;
    padding: 2% 0;
} 
.goback:hover + .instruction {
    display: flex;
    padding: 5px;
    width: 50px;
    font-size: 10px;
    text-align: center;
    color: rgba(60, 60, 60, 0.66);
    background-color: rgba(247, 247, 247, 0.90);
    box-shadow: 0px 0px 1px 0px var(--color-background); 
    position: absolute;
    left: 40px;
    top: 5px;
} }

</style>