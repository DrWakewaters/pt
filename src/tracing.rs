#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub enum Tracing {
    Direct,
    Light,
    Importance,
    Bidirectional,
}
