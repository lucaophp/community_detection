object Caminhos{
	def distance(graph:Graph[Int,Int],vertexInitial:Long):Graph[Int,Array[(Long,Double)]]={
		val initialGraph = graph.mapVertices((id,vertice)=>{
				if(id==vertexInitial) (0.0,Array()) else (Double.PositiveInfinity)
			})
	}
}