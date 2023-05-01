public class genlist<T>{
	public T[] data;
	public int size => data.Length;
	public T this[int i] => data[i];
	public genlist(){ data = new T[0]; }
	public void add(T item){
		T[] newdata = new T[size+1];
		System.Array.Copy(data,newdata,size);
		newdata[size]=item;
		data=newdata;
	}
	public void print(string s){
		System.Console.WriteLine(s);
		for(int i=0; i<size; i++){
			System.Console.WriteLine(data[i]);
		}
	}
}
